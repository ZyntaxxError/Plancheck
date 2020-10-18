////////////////////////////////////////////////////////////////////////////////
//  SRTcheck.cs
//
//  ESAPI v15.5 Script for simple plan parameter checks
//  
////////////////////////////////////////////////////////////////////////////////

using System;
using System.Collections.Generic;
using System.Text;
using System.Windows;
using System.Linq;
using VMS.TPS.Common.Model.API;
using VMS.TPS.Common.Model.Types;
using System.IO;
using System.Diagnostics;
using System.CodeDom.Compiler;


/* TODO: create plan category (enum) and make category specific checks for tbi, sbrt, electron, etc
 * TODO: IsPlanSRT fails if a revision made...
 * TODO Plan sum; foreach plan: check, easier to do in build and wpf...
 * TODO SBRT kolla body vs skin och om vakumkudde exists
 * */

namespace VMS.TPS
{
	class Script
	{
		public Script()
		{
		}
		public void Execute(ScriptContext context)
		{

			if (context.Patient == null || context.PlanSetup == null)
			{
				MessageBox.Show("Please select a patient and plan in active context window.");
			}
			else
			{
				Patient patient = context.Patient;
				Course course = context.Course;
				PlanSetup plan = context.PlanSetup;
				//enum planCategory{
				string courseIntent = context.Course.Intent;
				string courseId = context.Course.Id;

				string messageTitle = "Quick check on " + courseId + " " + plan.Id + "\n\n";
				//CheckPlanCategory(plan) +
				string message =
				CheckCourseIntent(courseIntent, plan) + "\n" +
				CheckPlanNamingConvention(plan) +
				CheckClinProt(plan) +
				CheckTargetVolumeID(plan, plan.StructureSet) +
				CheckCouchStructure(plan.StructureSet) +
				CheckForBolus(plan) +
				CheckSetupField(plan) + "\n" +
				//"Treatment fields \n" +
				//"ID \t Energy \t Tech. \t Drate \t MLC \n" +
				CheckFieldRules(plan) +
				CheckFieldNamingConvention(course, plan) +
				CheckStructureSet(plan) +
				DeltaShiftFromOrigin(plan);
                //SimpleCollisionCheck(plan);

                if (message.Length == 0)
                {
					message = "No comments...";
                }

				MessageBox.Show(message, messageTitle);
			}
		}

        private string DeltaShiftFromOrigin(PlanSetup plan)
        {
			VVector deltaShift = 10*DeltaShiftInmm(plan, plan.StructureSet.Image.UserOrigin, plan.Beams.First().IsocenterPosition);

			string delta = "\n\niso-Lat: \t" + deltaShift.x.ToString("0.0") + " mm\n";
			delta += "iso-Vrt: \t" + deltaShift.y.ToString("0.0") + " mm\n";
			delta += "iso-Lng: \t" + deltaShift.z.ToString("0.0") + " mm\n\n\n";
			return delta;
		}

        private VVector DeltaShiftInmm(PlanSetup plan, VVector dicomOriginalPosition, VVector dicomFinalPosition)
        {
			Image image = plan.StructureSet.Image;
			VVector eclipseOriginalPosition = image.DicomToUser(dicomOriginalPosition, plan);
			VVector eclipseFinalPosition = image.DicomToUser(dicomFinalPosition, plan);
			VVector deltaShift = eclipseOriginalPosition - eclipseFinalPosition;

			// Have to change sign of Vrt before returning due to difference in coordinate system in Eclipse and the machine
			deltaShift.y *= -1;

			return deltaShift;
		}

        private string CheckForBolus(PlanSetup plan)
        {
			string cResults = string.Empty;
			
			int noOfBoluses = 0;
			// iterate through all beams in plan, and then all boluses in beam
            foreach (var beam in plan.Beams)
            {
                foreach (var bolus in beam.Boluses)
                {
					noOfBoluses += beam.Boluses.Count();
					cResults += bolus.Id;
					cResults += bolus.MaterialCTValue.ToString("0.0");
                }
            }
			cResults += "\n* Check if bolus should be linked to all fields. Nr of boluses: " +  noOfBoluses.ToString() + "\n";
			return cResults;
        }


		// TMI: order plans from head to toe (by isopos in dicom, reverse), check ID numbering
		// 

        private string CheckStructureSet(PlanSetup plan)
        {
			string cResults = string.Empty;


			StructureSet ss = plan.StructureSet;
			Image image = ss.Image;
			Structure skin = ss.Structures.Where(s => s.Id == "skinn").SingleOrDefault();
			if (skin == null || skin.IsEmpty)
			{
				cResults = "* missing structure Skin \n";
			}
			return cResults;
			// Search for structures obligatory for an SBRT plan
		}




        /*
		private string CheckPlanCategory(PlanSetup plan)
		{
			string cResults = "";
			if (IsPlanSRT(plan))
			{
				cResults = "Running tests assuming this is a SRT plan";
			}
			return cResults;
		}

		enum PlanCat
		{
			SBRT,
			Electron,
			CRT,
			IMRT,
			TBI,
			TMI,
		Mamill,
			Unknown
		}

		PlanCat planCat = PlanCat.SBRT;
		*/
        public void SetPlanCategory(Course course, PlanSetup plan)
		{
			// first the easy categories



		}

		// ********* 	Kontroll om SRT-plan (enbart baserad på fraktionering, aning klent men bör fungera) *********

		public bool IsPlanSRT(PlanSetup plan)
		{
			bool fractionationSRT = ((plan.NumberOfFractions > 2 && plan.DosePerFraction.Dose >= 8.0 && plan.TotalDose.Dose >= 40.0) || (plan.NumberOfFractions >= 5 && plan.DosePerFraction.Dose >= 6 && plan.TotalDose.Dose >= 35.0));
			// missade 7Gy x 5
			return fractionationSRT;
		}
		/*public bool IsPlanElectron(PlanSetup plan)
		{
			bool fractionationSRS = ((plan.NumberOfFractions > 2 && plan.DosePerFraction.Dose >= 8.0 && plan.TotalDose.Dose >= 40.0) || (plan.NumberOfFractions > 7 && plan.DosePerFraction.Dose >= 6 && plan.TotalDose.Dose >= 45.0));

			return fractionationSRS;
		}*/




		// ********* Helper method for checking iso coordinates, returns VVector  
		// transforms coordinates from dicom to coordinates based on coronal view from table end, dicom-origo the same (usually center of image in Lat, below table in vrt)
		// TODO: check if this clashes with other properties in VVector   TODO: check if all fields same isocenter
		// TODO: would be nice if coordinates instead originates from center of image (or center of couch) in Lat, and Couch top surface in Vrt (this would also mean a 
		// chance to predict/estimate absolute couch coordinates in lat and vrt)


		// overloaded method , should keep only one of them and instead check if all beams have the same isocenter
		public VVector IsoPositionFromTableEnd(Beam beam, string treatOrient)
		{
			VVector beamIso = beam.IsocenterPosition; // mm from Dicom-origo
			switch (treatOrient)
			{
				case "FeetFirstSupine":
					{
						beamIso.x *= -1;
						break;
					}
				case "FeetFirstProne":
					{
						beamIso.y *= -1;
						break;
					}
				case "HeadFirstProne":
					{
						beamIso.x *= -1;
						beamIso.y *= -1;
						break;
					}
				default:
					break;
			}
			return beamIso;
		}



		// ********* Kontroll av Course Intent, kollar om ifyllt, annars enbart SRS/SBRT-planer  *********

		public string CheckCourseIntent(string cIntent, PlanSetup plan)
		{
			string cResults = "";
			if (cIntent.Trim().Length == 0)
			{
				cResults = "** Course intent is empty!";
			}
			else if (IsPlanSRT(plan) && !cIntent.Equals("SRT"))
			{
				cResults += "** Change course intent to SRT!";
			}
			return cResults + "\n";
		}

		// ********* 	Kontroll av kliniskt protokoll-ID om sådant kopplat till plan och där antal fraktioner ges av protokoll-ID, jämförs med plan-fraktionering *********
		// finns protokoll som inte har antal fraktioner i ID, går ej testa (finns metod, kolla detta) TODO

		public string CheckClinProt(PlanSetup plan)
		{
			string cResults = "";
			int fractionsInProtocol = 0;
			if (plan.ProtocolID.Length != 0)
			{
				int protocolFractionIndex = plan.ProtocolID.IndexOf('#'); // find the index of the symbol indicating nr of fractions
				if (protocolFractionIndex != -1)                            // if there are no fraction specification in ID skip the test
				{
					string protocolFrNrInfo = plan.ProtocolID.Substring(protocolFractionIndex - 2, 2).Trim();       // retrieve the two characters before the #, and remove whitespaces
					if (Int32.TryParse(protocolFrNrInfo, out fractionsInProtocol))                  // try parsing it to int and, if successful, compare to plan fractions
					{
						if (fractionsInProtocol != plan.NumberOfFractions)
						{
							cResults = "** Check the attached clinical protocol! \n \n";
						}
					}
				}

			}
			return cResults;
		}

		// ********* 	Targetvolym; kollar att det är valt och av Dicom-typen PTV *********

		public string CheckTargetVolumeID(PlanSetup plan, StructureSet sSet)
		{
			string cResults = "";
			if (string.IsNullOrEmpty(plan.TargetVolumeID))
			{
				cResults = "** No plan target volume selected \n";
			}
			else
			{
				// Search for structure in structure set with same id as target volume and checks if type is PTV, defaults to null if criteria not met
				Structure target = sSet.Structures.Where(s => s.Id == plan.TargetVolumeID).Where(s => s.DicomType == "PTV").SingleOrDefault();
				if (target == null)
				{
					cResults = "* Plan target volume should be of type PTV (ignore this if only field limits drawn) \n";
				}
				else
				{
					cResults = "Plan target volume: " + target.Id + "\n";
				}
			}
			return cResults + "\n";
		}



		// ********* 	Plannamn; kollar att det följer regler och att det inte är använt på behandlad eller godkänd plan i samma course *********
		// TODO: assumes no more than 9 plans in a single course, might break if TMI... in that case use iteration or RegEx
		private string CheckPlanNamingConvention(PlanSetup plan)
		{
			string cResults = string.Empty;
			
			char planIdFirstChar = plan.Id[0];                
			char planIdSecondChar = plan.Id[1];
			int planNumber;

			if (char.ToUpperInvariant(planIdFirstChar) == 'P'  && int.TryParse(planIdSecondChar.ToString(), out planNumber ))
            {

				CheckOtherPlanIdInCourse(plan, planNumber, ref cResults);

				if (!char.IsUpper(planIdFirstChar))
				{
					cResults += "* Please use upper case 'P' in plan ID.\n\n";
				}
			}
			return cResults;
		}

        private void CheckOtherPlanIdInCourse(PlanSetup plan, int planNumber, ref string cResults)
        {
			//TODO : make the function general so it can be used e.g. in field ID naming check
			//TODO: check for revision, automatic or manually made (i.e. former plan retired or completed early) and in that case the appended ':[revision number]'
			Course course = plan.Course;
			char planIdFirstChar = plan.Id[0];                
			char planIdSecondChar = plan.Id[1];
			int usedNumber;

			foreach (var ps in course.PlanSetups.Where(p => p.Id != plan.Id))
            {
				planIdFirstChar = ps.Id[0];
				planIdSecondChar = ps.Id[1];
				if (char.ToUpperInvariant(planIdFirstChar) == 'P' && int.TryParse(planIdSecondChar.ToString(), out usedNumber))
				{

					if (ps.ApprovalStatus == PlanSetupApprovalStatus.PlanningApproved || ps.ApprovalStatus == PlanSetupApprovalStatus.TreatmentApproved ||
					ps.ApprovalStatus == PlanSetupApprovalStatus.Completed || ps.ApprovalStatus == PlanSetupApprovalStatus.CompletedEarly)
                    {
                        if (planNumber == usedNumber)
                        {
							cResults += "* Plan ID '" + planIdFirstChar + planIdSecondChar + "' has already been used in a " + ps.ApprovalStatus + " plan:\n" + ps.Id + "\n\n";
						}
                    }
				}
			}
		}


        // ********* 	Kontroll av att bordsstruktur existerar, är av rätt typ, inte är tom och har korrekt HU 	********* 
        // begränsningar: kollar ej positionering i förhållande till body, ej heller att den inte är kapad. Rätt bordstyp kollas enbart på namn
        // Exact IGRT Couch, thin, medium och thick kan i princip vara korrekt beroende på lokalisation...

        public string CheckCouchStructure(StructureSet SSet)
		{
			string cResult = "";
			bool couchExt = false;
			bool couchInt = false;
			foreach (Structure s in SSet.Structures)
			{
				if (s.Id.Contains("CouchSurf") && !s.IsEmpty && s.Name.Contains("Exact IGRT Couch") && s.DicomType == "SUPPORT")
				{
					double couchExtHU;
					s.GetAssignedHU(out couchExtHU);
					if (Math.Round(couchExtHU) == -300)
					{
						couchExt = true;
					}
				}
				if (s.Id.Contains("CouchInt") && !s.IsEmpty && s.Name.Contains("Exact IGRT Couch") && s.DicomType == "SUPPORT")
				{
					double couchIntHU;
					s.GetAssignedHU(out couchIntHU);
					if (Math.Round(couchIntHU) == -1000)
					{
						couchInt = true;
					}
				}
			}
			if (!couchExt || !couchInt)
			{
				cResult = "** Check couch structure! \n";
			}
			return cResult;
		}


        #region Treatment field naming convention

        // ********* 	Kontroll att numrering av fält är konsekutivt, och att inte fältnumret använts i någon godkänd eller behandlad plan i samma Course *********
        // kollar dock ej om man hoppar över ett nummer mellan planer
        // TODO: recommend beam number change, prereq; need to check plan status and plan ID numbers
		// OK to to if plan numbers differ
        private string CheckFieldNamingConvention(Course course, PlanSetup plan)
		{
			string cResult = string.Empty;
			int tempNumber = 1000;
			int smallestBeamNumber = tempNumber;
			int number = tempNumber;
			var beamNumbersInPlan = new List<int>();

			// neccessary to check if Id is an integer first and find the smallest
			foreach (var beam in plan.Beams.Where(b => !b.IsSetupField).OrderBy(b => b.Id))
			{
				if (Int32.TryParse(beam.Id.Trim(), out number))
				{
					if (number < smallestBeamNumber)
					{
						smallestBeamNumber = number;
					}
					beamNumbersInPlan.Add(number);
				}
				else
				{
					cResult = " * Check field naming convention, field ID should be an integer  \t";
					beamNumbersInPlan.Clear();
					break;
				}
			}

			if (beamNumbersInPlan.Count != 0)
			{
				beamNumbersInPlan = beamNumbersInPlan.OrderBy(b => b).ToList();
				CheckConsecutiveFieldNumbers(beamNumbersInPlan, smallestBeamNumber, ref cResult);
				CheckIfBeamNumberUsedInOtherPlans(course, plan, beamNumbersInPlan, ref cResult);
			}

			return cResult;
		}


		//Check that the beam numbers within the plan are consecutive 
		private void CheckConsecutiveFieldNumbers(List<int> beamNumbersInPlan, int smallestBeamNumber, ref string cResult)
		{
			for (int i = 0; i < beamNumbersInPlan.Count; i++)
			{
				if (beamNumbersInPlan[i] == smallestBeamNumber + i)
				{
					i++;
				}
				else
				{
					cResult += " * Fields should be consecutively named \n";
					break;
				}
			}
		}

		// Check that beam number doesn't exist in other approved or completed plans within the same Course (Retired is ok if revision...)
		// TODO: fails to check if a number is skipped... however, that should only be done if the plan is approved 
		// TODO: check all plans regardless of status for plans where plan ID number is smaller than the plan under check
		private void CheckIfBeamNumberUsedInOtherPlans(Course course, PlanSetup plan, List<int> beamNumbersInPlan, ref string cResult)
		{
			List<int> beamNumbersInOtherPlans = new List<int>();
			List<int> beamNumbersUsed = new List<int>();
			foreach (var ps in course.PlanSetups.Where(p => p.Id != plan.Id))
            {
				beamNumbersInOtherPlans = GetBeamNumbersFromOtherPlans(ps);
				foreach (var n in beamNumbersInPlan)
                {
					if (beamNumbersInOtherPlans.Contains(n))
					{
						beamNumbersUsed.Add(n);
					}
				}
                if (beamNumbersUsed.Count != 0)
                {
					cResult += " * Field ID already used in a " + ps.ApprovalStatus + " plan: \t" + ps.Id + " ( ";
                    foreach (var fieldId in beamNumbersUsed)
                    {
						cResult += fieldId.ToString() + ", ";
					}
					cResult = cResult.Remove(cResult.Length - 2, 1);
					cResult += ")\n";
				}
				beamNumbersInOtherPlans.Clear();
				beamNumbersUsed.Clear();
			}
		}

		/// <summary>
		/// Get beam numbers from plans other than "plan" in the same course as "plan"
		///	TODO: excludes "retired", this might be ok if revision and plan name is the same but with appended revision number
		/// </summary>
		/// <param name="plan"></param>
		/// <returns></returns>

		private List<int> GetBeamNumbersFromOtherPlans(PlanSetup plan)
		{
			var beamNumbers = new List<int>();
			int number = 1000;

			if (plan.ApprovalStatus == PlanSetupApprovalStatus.PlanningApproved || plan.ApprovalStatus == PlanSetupApprovalStatus.TreatmentApproved || 
				plan.ApprovalStatus == PlanSetupApprovalStatus.Completed || plan.ApprovalStatus == PlanSetupApprovalStatus.CompletedEarly)
            {
				foreach (var beam in plan.Beams.Where(b => !b.IsSetupField).OrderBy(b => b.Id))
				{
					if (Int32.TryParse(beam.Id.Trim(), out number))
					{
						beamNumbers.Add(number);
					}
				}
            }
			return beamNumbers;
		}


		#endregion Treatment field naming convention


		// ********* 	Kontroll av Setup-fält; namngivning och ej bordsvridning ********* 

		public string CheckSetupField(PlanSetup plan)
		{
			string cResults = "";
			int countSetupfields = 0;
			foreach (var beam in plan.Beams)
			{
				if (beam.IsSetupField)
				{
					cResults = cResults + "Setup-field: \t" + beam.Id + "\t \t";
					countSetupfields++;
					if (beam.ControlPoints.First().PatientSupportAngle != 0)
					{
						cResults += "** Couch angle not 0!";
					}
					if (beam.Id.ToUpper().Substring(0, 2).Equals(plan.Id.ToUpper().Substring(0, 2))) // naming convention, should start with first two char in plan ID
					{
						if (beam.Id.ToUpper().Contains("CBCT")) // no extra checks for cbct-setup field neccessary
						{
							cResults = cResults + "OK" + "\n";
						}
						else
						{
							cResults = cResults + CheckPlanarSetupFields(beam, plan) + "\n";    // extra checks for planar setup fields
						}
					}
					else
					{
						cResults = cResults + "** ID should start with Plan-ID" + "\n";
					}
				}
			}
			if (countSetupfields == 0)
			{
				cResults = cResults + "missing!" + "\t \t" + "** Insert Setup field" + "\n";
			}
			return cResults;
		}

		// ********* 	Kontroll av planara Setup-fält (ej namngivna cbct), namngivning efter gantryvinkel	********* 

		public string CheckPlanarSetupFields(Beam beam, PlanSetup plan)
		{
			string cResults = "";
			string trimmedID = beam.Id.Substring(2).Trim();         // start iteration at index 2 (PX index 0 and 1)
			int gantryAngleInBeamID = 1000;                         // dummy number
			int latLimitForCollisionCheck = 40;
			int test = 1000;
			for (int i = 1; i < trimmedID.Length + 1; i++)
			{
				if (Int32.TryParse(trimmedID.Substring(0, i), out test))                    //  step up one char at a time and try parsing it to int. If successful assign it to gantryAngleInBeamID
				{
					gantryAngleInBeamID = test;
				}
			}
			if (gantryAngleInBeamID != 1000)
			{
				if (gantryAngleInBeamID == Math.Round(beam.ControlPoints.First().GantryAngle))
				{
					double lat = IsoPositionFromTableEnd(beam, plan.TreatmentOrientation.ToString()).x;
					if (Math.Abs(lat) < latLimitForCollisionCheck)
					{
						cResults += "OK";
					}
					else
					{
						// make a simple collision check based on isocenter coordinates
						bool tableToRight = (lat < 0);  // table move to right when viewed from table end
						switch (Convert.ToInt32(beam.ControlPoints.First().GantryAngle))
						{
							case 180:
								{
									if (tableToRight)
									{
										cResults += "OK";
									}
									else
									{
										cResults += "* Check for collision";
									}
									break;
								}
							case 90:
								{
									if (tableToRight)
									{
										cResults += "* Check for collision";
									}
									else
									{
										cResults += "OK";
									}
									break;
								}
							case 0:
								{
									if (tableToRight)
									{
										cResults += "* Check for collision";
									}
									else
									{
										cResults += "OK";
									}
									break;
								}
							default:
								cResults += "OK";
								break;
						}
					}
				}
				else
				{
					cResults += "* Check name!";
				}
			}
			else
			{
				cResults += "* Check name!";
			}
			return cResults;
		}


		// ********* 	Kontroll av diverse fältregler och "best practices"	********* 
		// TODO: better sorting, maybe one general and then divided in categories (enum). Missing case for Static-I and MLC doseDynamic (IMRT)

		public string CheckFieldRules(PlanSetup plan)
		{
			string cResults = "";
			string remarks = "";
			int countFields = plan.Beams.Count();
			int countSetupFields = plan.Beams.Where(b => b.IsSetupField).Count();
			int countTreatFields = countFields - countSetupFields;
			int countSRSRemarks = 0;                    // All fields should be SRS if SBRT- or SRS-plan
			int countDoseRateFFFRemarks = 0;            // Dose rate should be maximum for FFF
			int countArcDynFFFRemarks = 0;              // The energy should be FFF if dynamic arc used and SRS
			int countMLCStaticFFFRemarks = 0;           // The energy should be FFF if static MLC (regardless of arc or static field) used and SRS
			int countArcDynCollAngleRemarks = 0;        // The collimator angle should be between +/-5 deg if dynamic arc used
			int countArcCW = 0;
			int countArcCCW = 0;                        // the absolute difference between CW and CCW should be less than two...

			foreach (var beam in plan.Beams.Where(b => !b.IsSetupField).OrderBy(b => b.Id))
			{
				cResults = cResults + beam.Id + "\t" + beam.EnergyModeDisplayName + "\t" + beam.Technique.Id + "\t" + beam.DoseRate + "\t" + beam.MLCPlanType + "\n";
				if (!beam.Technique.Id.Contains("SRS") && IsPlanSRT(plan))
				{
					if (countSRSRemarks < 1) { remarks = remarks + "** Change technique to SRS-" + beam.Technique.Id + "! \n"; };
					countSRSRemarks++;
				}
				if (beam.EnergyModeDisplayName.Contains("FFF"))
				{
					remarks += CheckDoseRateFFF(beam, ref countDoseRateFFFRemarks);
				}
				if (beam.MLCPlanType == MLCPlanType.ArcDynamic && IsPlanSRT(plan))
				{
					remarks += CheckArcDynFFF(beam, ref countArcDynFFFRemarks);
					remarks += CheckArcDynCollAngle(beam, ref countArcDynCollAngleRemarks);
				}
				if (beam.MLCPlanType == MLCPlanType.Static && IsPlanSRT(plan))
				{
					remarks += CheckMLCStaticFFF(plan, beam, ref countMLCStaticFFFRemarks);
				}
				// the absolute difference between CW and CCW should be less than two...
				if (beam.GantryDirection == GantryDirection.CounterClockwise)
				{
					countArcCCW++;
				}
				if (beam.GantryDirection == GantryDirection.Clockwise)
				{
					countArcCW++;
				}
				if (countTreatFields == 1 && beam.Technique.Id.Contains("ARC"))
				{
					remarks += CheckArcStartStop(beam);
				}
			}
			if (Math.Abs(countArcCCW - countArcCW) > 1)
			{
				remarks += "** Check the arc directions! \t";
			}
			return cResults + "\n" + remarks;
		}

		// ********* 	Enkel kontroll av eventuell kollisionsrisk, enbart planara setupfält	*********
		/*private string SimpleCollisionCheck(PlanSetup plan)
		{
			string cResults = "";
            if (IsoPositionFromTableEnd(plan).x < -40 && )
            {

            }


			return cResults;
			// TODO: check if there are "extended" accessible or perhaps gantry direction even if its static
		}*/


		// ********* 	Kontroll av dosrat vid FFF	*********

		public string CheckDoseRateFFF(Beam beam, ref int countDoseRateFFFRemarks)
		{
			string cResults = "";

			if ((beam.DoseRate < 1400 && beam.EnergyModeDisplayName.Contains("6")) || (beam.DoseRate < 2400 && beam.EnergyModeDisplayName.Contains("10")))
			{
				if (countDoseRateFFFRemarks < 1)
				{
					cResults = "* Consider changing dose rate to maximum!\n";
					countDoseRateFFFRemarks++;
				}
			}
			return cResults;
		}

		// ********* 	Kontroll av energi vid Dynamic Arc och SRT	*********

		public string CheckArcDynFFF(Beam beam, ref int countArcDynFFFRemarks)
		{
			string cResults = "";
			if (countArcDynFFFRemarks < 1 && !beam.EnergyModeDisplayName.Contains("FFF"))
			{
				cResults = "* Consider changing energy to FFF \n";
				countArcDynFFFRemarks++;
			}
			return cResults;
		}

		// ********* 	Kontroll av energi vid Statisk MLC och SRT	*********

		public string CheckMLCStaticFFF(PlanSetup plan, Beam beam, ref int countMLCStaticFFFRemarks)
		{
			string cResults = "";
			// only if no wedges in plan
			if (countMLCStaticFFFRemarks < 1 && !beam.EnergyModeDisplayName.Contains("FFF") && !anyWedgesInPlan(plan))
			{
				cResults = "* Consider changing energy to FFF \n";
				countMLCStaticFFFRemarks++;
			}
			return cResults;
		}
		public bool anyWedgesInPlan(PlanSetup plan)
		{
			bool anyWedgesInPlan = false;
			foreach (var beam in plan.Beams.Where(b => !b.IsSetupField))
			{
				if (beam.Wedges.Any())
				{
					anyWedgesInPlan = true;
				}
			}
			return anyWedgesInPlan;
		}


		// ********* 	Kontroll av kollimatorvinkel vid Dynamic Arc	*********

		public string CheckArcDynCollAngle(Beam beam, ref int countArcDynCollAngleRemarks)
		{
			string cResults = "";
			if (countArcDynCollAngleRemarks < 1 && beam.ControlPoints.First().CollimatorAngle > 5.0 && beam.ControlPoints.First().CollimatorAngle < 355.0)
			{
				cResults = "** Collimator angle for DynArc should be between +/- 5 deg \n";
				countArcDynCollAngleRemarks++;
			}
			return cResults;
		}

		// ********* 	Kontroll av start och stop-vinkel vid Arc och enbart ett fält, bör helst (om inte fullvarv) stanna ovan pat ********

		public string CheckArcStartStop(Beam beam)
		{
			string cResults = "";
			double start = beam.ControlPoints.First().GantryAngle;
			double stop = beam.ControlPoints.Last().GantryAngle;
			if ((start > 270.0 || start < 90.0) && (stop > 90.0 && stop < 270.0))
			{
				cResults = "* Consider reversing arc direction to make the gantry stop above the patient when treatment is finished. \n";
			}
			return cResults;
		}





	}
}
