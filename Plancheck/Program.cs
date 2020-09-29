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


//using System.Windows.Forms;
//using System.Windows.Forms;



/* TODO: create plan category (enum) and make category specific checks for tbi, sbrt, electron, etc
 * TODO: IsPlanSRT fails if a revision made...
 * TODO Plan sum; foreach plan: check, easier to do in build and wpf...
 * 
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

				string message =
				"Quick check on " + courseId + " " + plan.Id + "\n\n" +
				CheckPlanCategory(plan) +
				CheckCourseIntent(courseIntent, plan) + "\n" +
				CheckClinProt(plan) +
				CheckPlanProp(plan, plan.StructureSet) +
				CheckCouchStructure(plan.StructureSet) +
				CheckSetupField(plan) + "\n" +
				"Treatment fields \n" +
				"ID \t Energy \t Tech. \t Drate \t MLC \n" +
				CheckFieldRules(plan) +
				CheckConsecutiveFieldNaming(course, plan);
				//SimpleCollisionCheck(plan);
				MessageBox.Show(message);

				if (IsPlanSRT(plan))
				{
					string checkSBRTOrigo = CheckUserOriginInSBRTFrame(plan);
					MessageBox.Show(checkSBRTOrigo);
					string checkSBRTiso = GetIsoCoordInSBRTFrame(plan);
					MessageBox.Show(checkSBRTiso);

					//string isoCoordInSBRT = GetIsoCoordInSBRTFrame(plan);
					// TODO check if isocenter in same plane as user origo, not neccesary though as there can be multiple isocenters (muliple plans) and there is no strict rule...
				}
			}
		}



		private string CheckPlanCategory(PlanSetup plan)
		{
			string cResults = "";
			if (IsPlanSRT(plan))
			{
				cResults = "Running tests assuming this is a SRT plan";
			}
			return cResults;
		}




		// ********* 	Kontroll om SRT-plan (enbart baserad på fraktionering, aning klent men bör fungera) *********

		public bool IsPlanSRT(PlanSetup plan)
		{
			bool fractionationSRT = ((plan.NumberOfFractions > 2 && plan.DosePerFraction.Dose >= 8.0 && plan.TotalDose.Dose >= 40.0) || (plan.NumberOfFractions > 7 && plan.DosePerFraction.Dose >= 6 && plan.TotalDose.Dose >= 45.0));

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
		public VVector IsoPositionFromTableEnd(PlanSetup plan)
		{
			//var image = plan.StructureSet.Image;
			VVector planIso = plan.Beams.First().IsocenterPosition; // mm from Dicom-origo
			switch (plan.TreatmentOrientation.ToString())
			{
				case "FeetFirstSupine":
					{
						planIso.x *= -1;
						break;
					}
				case "FeetFirstProne":
					{
						planIso.y *= -1;
						break;
					}
				case "HeadFirstProne":
					{
						planIso.x *= -1;
						planIso.y *= -1;
						break;
					}
				default:
					break;
			}
			return planIso;
		}

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
				if (protocolFractionIndex != -1)                            // if there are no fraktion specification in ID skip the test
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

		public string CheckPlanProp(PlanSetup plan, StructureSet sSet)
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

		// ********* 	Kontroll att numrering av fält är konsekutivt, och att inte fältnumret använts i någon godkänd eller behandlad plan i samma Course *********
		// kollar dock ej om man hoppar över ett nummer mellan planer

		private string CheckConsecutiveFieldNaming(Course course, PlanSetup plan)
		{
			string cResult = "";

			// ugly code, but works, refactor somehow... but neccessary to check if Id is a number first and find the smallest

			int smallestBeamNumber = 1000;
			int number = 1000;
			var beamNumbersInPlan = new List<int>();
			var beamNumbersInOtherPlans = new List<int>();
			beamNumbersInOtherPlans = GetBeamNumbersFromOtherPlans(course, plan);
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
					cResult = " * Check field naming convention \t";
					beamNumbersInPlan.Clear();
					break;
				}
			}
			//Check that the beam numbers within the plan are consecutive and doesn't exist in other approved or completed plans within the same Course (Retired is ok if revision...)
			// hmm, failes to check if a number is skipped... however, that should only be done if the plan is approved TODO
			int i = 0;
			if (smallestBeamNumber != 1000)
			{
				foreach (var n in beamNumbersInPlan.OrderBy(b => b))
				{
					if (n == smallestBeamNumber + i)
					{
						if (beamNumbersInOtherPlans.Contains(n))
						{
							cResult = " * Beam number used in other plan  \t";
							break;
						}
						i++;
					}
					else
					{
						cResult = " * Fields should be consecutively named \t";
						break;
					}
				}
			}
			return cResult;
		}

		/// <summary>
		/// Get beam numbers from plans other than "plan" in the same course as "plan"
		/// </summary>
		/// <param name="course"></param>
		/// <param name="plan"></param>
		/// <returns></returns>
		private static List<int> GetBeamNumbersFromOtherPlans(Course course, PlanSetup plan)
		{
			var beamNumbers = new List<int>();
			int number = 1000;
			foreach (var ps in course.PlanSetups.Where(p => p.Id != plan.Id))
			{
				if (ps.ApprovalStatus.ToString().Equals("PlanningApproved") || ps.ApprovalStatus.ToString().Equals("TreatmentApproved") || ps.ApprovalStatus.ToString().Equals("Completed") || ps.ApprovalStatus.ToString().Equals("CompletedEarly"))
				{
					foreach (var beam in ps.Beams.Where(b => !b.IsSetupField).OrderBy(b => b.Id))
					{
						if (Int32.TryParse(beam.Id.Trim(), out number))
						{
							beamNumbers.Add(number);
						}
					}
				}
			}
			return beamNumbers;
		}

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

		//*********HELPER METHODS**************

		bool SameSign(double num1, double num2)
		{
			return num1 >= 0 && num2 >= 0 || num1 < 0 && num2 < 0;
		}

		public class PatternGradient
		{
			public List<double> DistanceInMm { get; set; }
			public List<int> GradientHUPerMm { get; set; }
			public List<int> PositionToleranceMm { get; set; }
			public int GradIndexForCoord { get; set; }
		}

		/// <summary>
		/// 
		/// </summary>
		/// <param name="plan"></param>
		/// <returns></returns>
		public string CheckUserOriginInSBRTFrame(PlanSetup plan)
		{
			string userOrigoCheck;
			int userOrigoLongSRS = 0;
			var image = plan.StructureSet.Image;
			double bottom;
			double left;				// left and right side only used to get longcoordinates
			double right;

			var coord = GetTransverseCoordInSBRTFrame(plan, image.UserOrigin);
			left = coord[0];
			right = coord[1];
			bottom = coord[2];
			double lateralCenterSBRT = (left+right)/2;
			{
				double coordSRSLong = GetSRSLongCoord(plan, image.UserOrigin, (int)Math.Round(bottom), left, right);
				if (coordSRSLong != 0)
				{
					userOrigoLongSRS = (int)Math.Round(coordSRSLong);
				}
			}
			int userOrigoVrtSRS = Convert.ToInt32(Math.Round(Math.Abs(image.UserOrigin.y - bottom)));
			int userOrigoLatSRS = Convert.ToInt32(Math.Round(-(image.UserOrigin.x - lateralCenterSBRT - 300)));		// TODO: this works for HFS and FFS. HFP and FFP unhandled
			
			if (bottom == 0 || userOrigoLongSRS == 0)
			{
				userOrigoCheck = "Cannot find the SRS-frame, no automatic check of User origo possible.";
			}
			else if (Math.Abs(image.UserOrigin.y - (bottom - 95)) < 3 && Math.Abs(image.UserOrigin.x - lateralCenterSBRT) < 3)
			{
				userOrigoCheck = "Estimated position of user origin in SRS frame coordinates from image profiles: \n\n" +
				" Lat: " + userOrigoLatSRS + "\t Vrt: " + userOrigoVrtSRS + "\t Lng: " + userOrigoLongSRS + "\t (+/- 2mm)" +
				"\n\n";
			}
			else
			{
				userOrigoCheck = "Check position of user origo";
			}
			return userOrigoCheck;
		}


		public string GetIsoCoordInSBRTFrame(PlanSetup plan)
		{
			string isoSBRTresults;
			int isoLongSRS = 0;
			var image = plan.StructureSet.Image;
			double bottom;
			double left;                // left and right side only used to get longcoordinates
			double right;

			var coord = GetTransverseCoordInSBRTFrame(plan, plan.Beams.First().IsocenterPosition);
			left = coord[0];
			right = coord[1];
			bottom = coord[2];
			double lateralCenterSBRT = (left + right) / 2;
			{
				double coordSRSLong = GetSRSLongCoord(plan, plan.Beams.First().IsocenterPosition, (int)Math.Round(bottom), left, right);
				if (coordSRSLong != 0)
				{
					isoLongSRS = (int)Math.Round(coordSRSLong);
				}
			}
			int isoVrtSBRT = Convert.ToInt32(Math.Round(Math.Abs(plan.Beams.First().IsocenterPosition.y - bottom)));
			int isoLatSBRT = Convert.ToInt32(Math.Round(-(plan.Beams.First().IsocenterPosition.x - lateralCenterSBRT - 300)));     // TODO: this works for HFS and FFS. HFP and FFP unhandled

			// Get Vrt and Lat calculated directly from user origo for comparison
			int isoVrtSBRTFromUO = (int)(Math.Round(Math.Abs(plan.Beams.First().IsocenterPosition.y - image.UserOrigin.y - 95)));	// TODO: this works for HFS and FFS. HFP and FFP unhandled
			int isoLatSBRTFromUO = (int)(Math.Round(-(plan.Beams.First().IsocenterPosition.x - image.UserOrigin.x - 300)));         // TODO: this works for HFS and FFS. HFP and FFP unhandled

			if (bottom == 0 && isoLongSRS == 0)
			{
				isoSBRTresults = "Cannot find the SBRT frame, no automatic check of isocenter possible. Lateral and vertical calculated from user origo in paranthesis, assuming user origo is correctly positioned in Lat 300 and Vrt 95: \n\n" +
									"\n (Lat: " + isoLatSBRTFromUO + ")\t (Vrt: " + isoVrtSBRTFromUO + ")";
			}
            else if (isoLongSRS == 0)
            {
				isoSBRTresults = "Cannot find the SBRT frame long coordinate with profiles. Estimated position of isocenter in SBRT frame coordinates in transverse direction from image profiles (calculated from user origo in paranthesis): \n\n" +
				" Lat: " + isoLatSBRT + "\t Vrt: " + isoVrtSBRT + "\t (+/- 3mm)" +
				"\n (Lat: " + isoLatSBRTFromUO + ")\t (Vrt: " + isoVrtSBRTFromUO + ")";
			} 
			else
			{
				isoSBRTresults = "Estimated position of isocenter in SBRT frame coordinates from image profiles (calculated from user origo in paranthesis): \n\n" +
				" Lat: " + isoLatSBRT + "\t Vrt: " + isoVrtSBRT + "\t Lng: " + isoLongSRS + "\t (+/- 3mm)" +
				"\n (Lat: " + isoLatSBRTFromUO + ")\t (Vrt: " + isoVrtSBRTFromUO + ")";
			}
			return isoSBRTresults;
		}










		/// <summary>
		/// Gets frame of reference in the SBRT frame by returning the dicom coordinate of the SBRT frame bottom (represents 0 in frame coordinates)
		/// and left and right side of the frame in lateral direction, in the plane given by "dicomPosition"
		/// </summary>
		/// <param name="plan"></param>
		/// <param name="dicomPosition"></param>
		/// <returns></returns>
		public double[] GetTransverseCoordInSBRTFrame(PlanSetup plan, VVector dicomPosition)
		{
			VVector frameOfRefSBRT = dicomPosition;

			// get the dicom-position representing vertical coordinate 0 in the SBRT frame
			frameOfRefSBRT.y = GetSBRTBottomCoord(plan, dicomPosition);    // igores the position of dicomPosition in x and y and takes only z-position, takes bottom and center of image
			VVector frameOfRefSBRTLeft = frameOfRefSBRT;					// TODO: left and right designation depending of HFS FFS etc...
			VVector frameOfRefSBRTRight = frameOfRefSBRT;
			double[] returnCoo = new double[3];
			var image = plan.StructureSet.Image;

			// If bottom found, dubblecheck that the frameOfRefSBRT.y really is the bottom of the SBRT-frame by taking profiles in the sloping part of the SBRT frame and comparing
			// width with expected width at respective height above the bottom
			if (frameOfRefSBRT.y != 0 && DoubleCheckSBRTVrt(image, frameOfRefSBRT))
			{




                // debug 

				//int userOrigoVrtSRS = Convert.ToInt32(Math.Round(Math.Abs(image.UserOrigin.y - frameOfRefSBRT.y)));
				//MessageBox.Show("Found SBRT bottom  " + userOrigoVrtSRS);

				
				double[] coordSRSLat = new double[2];
				coordSRSLat = GetSRSLatCoord(plan, dicomPosition, (int)Math.Round(frameOfRefSBRT.y)); 
				frameOfRefSBRTLeft.x = coordSRSLat[0];
				frameOfRefSBRTRight.x = coordSRSLat[1];
				frameOfRefSBRT.x = (frameOfRefSBRTLeft.x + frameOfRefSBRTRight.x) / 2;  // middle of the SBRT frame in Lat
				// what are the chances that the dicom coord in x actually is 0.0? probably small but need to handle this by setting vrt to 0, TODO: could perhaps use nullable ref types 
				// rough check of the found width of the SBRT frame, allow for some flex of the frame. hardcoded magic number... TODO: this will brake if the gradient number changes
				if (frameOfRefSBRTLeft.x == 0 || frameOfRefSBRTRight.x == 0 || Math.Abs(frameOfRefSBRTLeft.x - frameOfRefSBRTRight.x) < 439 || Math.Abs(frameOfRefSBRTLeft.x - frameOfRefSBRTRight.x) > 445)
				{
					frameOfRefSBRT.x = 0;
					frameOfRefSBRT.y = 0;
				}
				MessageBox.Show("Frame width:  " + Math.Abs(frameOfRefSBRTLeft.x - frameOfRefSBRTRight.x).ToString("0.0"));
			}
			else
			{
				MessageBox.Show("Could not find the SBRT frame bottom. ");
			}

			returnCoo[0] = frameOfRefSBRTLeft.x;
			returnCoo[1] = frameOfRefSBRTRight.x;
			returnCoo[2] = frameOfRefSBRT.y;
			return returnCoo;

			// TODO check if isocenter in same plane as user origo, not neccesary though as there can be multiple isocenters (muliple plans) and there is no strict rule...
		}

        private bool DoubleCheckSBRTVrt(Image image, VVector frameOfRefSBRT)
        {
			// 10 mm above the bottom the expected width of the frame is 349 mm (third gradient, i.e. the inner surface of the inner wall) Changes fast with height... approx. 140 deg angle of box, 12 mm above; 356...
			bool checkResult = false;
			double xLeftUpperCorner = image.Origin.x - image.XRes / 2;  // Dicomcoord in upper left corner ( NOT middle of voxel in upper left corner)
			VVector leftProfileStart = frameOfRefSBRT;                // only to get the z-coord of the passed in VVector, x and y coord will be reassigned
			VVector rightProfileStart = frameOfRefSBRT;               // only to get the z-coord of the passed in VVector, x and y coord will be reassigned
			leftProfileStart.x = xLeftUpperCorner + image.XRes;         // start 1 pixel in left side
			rightProfileStart.x = xLeftUpperCorner + image.XSize * image.XRes - image.XRes;         // start 1 pixel in right side
			leftProfileStart.y = frameOfRefSBRT.y - 10;                 // 10 mm from assumed bottom    
			rightProfileStart.y = leftProfileStart.y;
			double stepsX = image.XRes;             //   (mm/voxel) to make the steps 1 pixel wide, can skip this if 1 mm steps is wanted

			VVector leftProfileEnd = leftProfileStart;
			VVector rightProfileEnd = rightProfileStart;
			leftProfileEnd.x += 200 * stepsX;                   
			rightProfileEnd.x -= 200 * stepsX;

			var samplesX = (int)Math.Ceiling((leftProfileStart - leftProfileEnd).Length / stepsX);

			var profLeft = image.GetImageProfile(leftProfileStart, leftProfileEnd, new double[samplesX]);
			var profRight = image.GetImageProfile(rightProfileStart, rightProfileEnd, new double[samplesX]);

			List<double> valHULeft = new List<double>();
			List<double> cooLeft = new List<double>();
			for (int i = 0; i < samplesX; i++)
			{
				valHULeft.Add(profLeft[i].Value);
				cooLeft.Add(profLeft[i].Position.x);
			}

			List<double> valHURight = new List<double>();
			List<double> cooRight = new List<double>();
			for (int i = 0; i < samplesX; i++)
			{
				valHURight.Add(profRight[i].Value);
				cooRight.Add(profRight[i].Position.x);
			}


			//***********  Gradient patter describing expected profile in HU of the Lax-box slanted side, from outside to inside **********

			PatternGradient sbrtSide = new PatternGradient();
			sbrtSide.DistanceInMm = new List<double>() { 0, 2, 20 };     // distance between gradients, mean values from profiling 10 pat 
			sbrtSide.GradientHUPerMm = new List<int>() { 100, -100, 100 };
			sbrtSide.PositionToleranceMm = new List<int>() { 0, 3, 3 };                        // tolerance for the gradient position
			sbrtSide.GradIndexForCoord = 2;                      // index of gradient position to return (zero based index), i.e. the start of the inner wall

			double[] coordBoxLat = new double[2];
			coordBoxLat[0] = GetCoordinates(cooLeft, valHULeft, sbrtSide.GradientHUPerMm, sbrtSide.DistanceInMm, sbrtSide.PositionToleranceMm, sbrtSide.GradIndexForCoord);
			coordBoxLat[1] = GetCoordinates(cooRight, valHURight, sbrtSide.GradientHUPerMm, sbrtSide.DistanceInMm, sbrtSide.PositionToleranceMm, sbrtSide.GradIndexForCoord);
			//coordBoxLat[2] = ((coordBoxRight + coordBoxLeft) / 2);
			if (coordBoxLat[0] != 0 && coordBoxLat[1] != 0 && Math.Abs(coordBoxLat[1]- coordBoxLat[0]) < 356 && Math.Abs(coordBoxLat[1] - coordBoxLat[0]) > 342)
			{
				MessageBox.Show("Width of slanted box side " + Math.Abs(coordBoxLat[1] - coordBoxLat[0]).ToString("0.0"));
				checkResult = true;
			}

			return checkResult;
        }










        /// <summary>
        /// Gets the coordinates of the bottom of the SBRT frame, given the plan and the position of interest
        /// Takes the position in center of image in z-coord given by "dicomPosition"
        /// </summary>
        /// <param name="plan"></param>
        /// <param name="dicomPosition"></param>
        /// <returns></returns>
        private int GetSBRTBottomCoord(PlanSetup plan, VVector dicomPosition)
		{
			//var planIso = plan.Beams.First().IsocenterPosition; // mm from Dicom-origo
			var image = plan.StructureSet.Image;
			//var imageUserOrigo = image.UserOrigin;             // mm från origo satt från CT vilket är dicom-origo!The user origin in DICOM coordinates in millimeter. 
			//var imageCTO = image.Origin;                // Ursprungligt origo satt från CT, mm från CT-origo TILL övre vänstra hörnet i första bilden!
			//The origin of the image. In other words, the DICOM coordinates of the center point of the upper-left hand corner voxel of the first image plane

			double imageSizeX = image.XRes * image.XSize;
			double imageSizeY = image.YRes * image.YSize;
			//double dist = VVector.Distance(imageCTO, imageUserOrigo);  // 3D distance from one koord to another, double 
			//var userIsoCoord = image.DicomToUser(planIso, plan);  //   coordinates from USER ORIGO to iso
			// VVector @ image center in iso plane in dicomcoordinates   What...?
			//VVector isoPlaneImageCenter = planIso;

			double xLeftUpperCorner = image.Origin.x - image.XRes / 2;  // Dicomcoord in upper left corner ( NOT middle of voxel in upper left corner)
			double yLeftUpperCorner = image.Origin.y - image.YRes / 2;  // Dicomcoord in upper left corner ( NOT middle of voxel in upper left corner)

			// instead of absolute image center => coord of first and last image voxel in x and y, XSize in voxels and XRes in mm/voxel
			// this seems like a smart idea if I can understand it...
			/*double xVoxStart = image.Origin.x;
			double xVoxEnd = image.Origin.x + image.XRes * image.XSize - image.XRes;
			double yVoxStart = image.Origin.y;
			double yVoxEnd = image.Origin.y + image.YRes * image.YSize - image.YRes;
			*/
			//**********  Image profile, need to define startpoint(VVector), endpoint(VVector) and number of samples (double[]) in that distance *********** 
			//In case of SBRT we want to get a profile in User origo or isocenter plane, middle of the image in x-direction and starting from bottom of the image in y.

			VVector bottomProfileStart = dicomPosition;                      // only to get the z-coord of the user origo, x and y coord will be reassigned
			bottomProfileStart.x = xLeftUpperCorner + imageSizeX / 2;           // center of the image in x-direction
			bottomProfileStart.y = yLeftUpperCorner + imageSizeY - image.YRes;  // start 1 pixel in from bottom...
			double steps = plan.StructureSet.Image.YRes;                        //  (mm/voxel) to make the steps 1 pixel wide, can skip this if 1 mm steps is wanted

			VVector bottomProfileEnd = bottomProfileStart;
			bottomProfileEnd.y -= 200 * steps;                                  // endpoint 200 steps in -y direction, i.e. 20 cm if 1 mm pixels

			var samplesY = (int)Math.Ceiling((bottomProfileStart - bottomProfileEnd).Length / steps);
			var profY = image.GetImageProfile(bottomProfileStart, bottomProfileEnd, new double[samplesY]);

			//***********  Gradient patter describing expected profile in HU of the sbrt-box bottom **********
			/*
			PatternGradient sbrt = new PatternGradient();
			sbrt.DistanceInMm = new List<double>() { 0, 4.4, 12.3 };//, 2 };        // distance between gradients, mean values from profiling 10 pat
			sbrt.GradientHUPerMm = new List<int>() { 100, -100, 100 };//, -100 };      // ,  inner shell can be separated from the box, larger tolerance
			sbrt.PositionToleranceMm = new List<int>() { 0, 2, 2 };//, 3 };                                   // tolerance for the gradient position
			sbrt.GradIndexForCoord = 2;   
			*/// index of gradient position to return (zero based index)

			PatternGradient sbrt = new PatternGradient();
			sbrt.DistanceInMm = new List<double>() { 0, 4.4, 12.3, 2 };        // distance between gradients, mean values from profiling 10 pat
			sbrt.GradientHUPerMm = new List<int>() { 90, -90, 90, -90 };      // ,  inner shell can be separated from the box, larger tolerance
			sbrt.PositionToleranceMm = new List<int>() { 0, 2, 1, 3 };        // tolerance for the gradient position, needs to be tight as some ct couch tops have almost the same dimensions
			sbrt.GradIndexForCoord = 2;                                     // index of gradient position to return (zero based index)





			// Imageprofile gets a VVector back, take the coordinates and respective HU and put them in two Lists of double, might be better ways of doing this...
			List<double> valHU = new List<double>();
			List<double> coo = new List<double>();
			string debug1 = "";
			for (int i = 0; i < samplesY; i++)
			{
				valHU.Add(profY[i].Value);
				coo.Add(profY[i].Position.y);
                if (profY[i].Value > -1000 && profY[i].Position.y <0)
                {
					debug1 += coo[i].ToString("0.0") + "\t" + (valHU[i]- valHU[i-1]).ToString("0.0") + "\n";
				}
			}
			//MessageBox.Show(debug1);

			// Get the coordinate (dicom) that represents inner bottom of SBRT frame 
			double coordBoxBottom = GetCoordinates(coo, valHU, sbrt.GradientHUPerMm, sbrt.DistanceInMm, sbrt.PositionToleranceMm, sbrt.GradIndexForCoord);
            // in the SBRT frame; VRT 0, which we are looking for, is approximately 1 mm above this gradient position, add 1 mm before returning
            if (coordBoxBottom != 0)
            {
				coordBoxBottom -= 1;        // TODO: this works for HFS and FFS, HFP and FFP should be handled
			}
			
			return (int)Math.Round(coordBoxBottom);
		}


		private double[] GetSRSLatCoord(PlanSetup plan, VVector dicomCoord, int coordSRSBottom)
		{

			// ************************************ get profiles in x direction, left and right side and determine center of box ********************

			var image = plan.StructureSet.Image;
			double xLeftUpperCorner = image.Origin.x - image.XRes / 2;  // Dicomcoord in upper left corner ( NOT middle of voxel in upper left corner)
			VVector leftProfileStart = dicomCoord;                // only to get the z-coord of the passed in VVector, x and y coord will be reassigned
			VVector rightProfileStart = dicomCoord;               // only to get the z-coord of the passed in VVector, x and y coord will be reassigned
			leftProfileStart.x = xLeftUpperCorner + image.XRes;         // start 1 pixel in left side
			rightProfileStart.x = xLeftUpperCorner + image.XSize * image.XRes - image.XRes;         // start 1 pixel in right side
			leftProfileStart.y = coordSRSBottom - 91.5;                 // hopefully between fidusles...     
			rightProfileStart.y = leftProfileStart.y;
			double stepsX = image.XRes;             //   (mm/voxel) to make the steps 1 pixel wide, can skip this if 1 mm steps is wanted

			VVector leftProfileEnd = leftProfileStart;
			VVector rightProfileEnd = rightProfileStart;
			leftProfileEnd.x += 100 * stepsX;                   // endpoint 100 steps in  direction
			rightProfileEnd.x -= 100 * stepsX;

			var samplesX = (int)Math.Ceiling((leftProfileStart - leftProfileEnd).Length / stepsX);

			var profLeft = image.GetImageProfile(leftProfileStart, leftProfileEnd, new double[samplesX]);
			var profRight = image.GetImageProfile(rightProfileStart, rightProfileEnd, new double[samplesX]);

			List<double> valHULeft = new List<double>();
			List<double> cooLeft = new List<double>();
			string debugLeft = "";
			for (int i = 0; i < samplesX; i++)
			{
				valHULeft.Add(profLeft[i].Value);
				cooLeft.Add(profLeft[i].Position.x);
				if (i > 0)
				{
					debugLeft += cooLeft[i].ToString("0.0") + "\t" + (valHULeft[i] - valHULeft[i - 1]).ToString("0.0") + "\n";
				}
			}


			List<double> valHURight = new List<double>();
			List<double> cooRight = new List<double>();
			//string debugRight = "";
			for (int i = 0; i < samplesX; i++)
			{
				valHURight.Add(profRight[i].Value);
				cooRight.Add(profRight[i].Position.x);
				//debugLeft += cooRight[i].ToString("0.0") + "\t" + valHURight[i].ToString("0.0") + "\n";
			}



			//MessageBox.Show(debugLeft);

			//***********  Gradient patter describing expected profile in HU of the Lax-box side, from outside to inside **********

			PatternGradient sbrtSide = new PatternGradient();
			sbrtSide.DistanceInMm = new List<double>() { 0, 2, 13 };     // distance between gradients, mean values from profiling 10 pat 
			sbrtSide.GradientHUPerMm = new List<int>() { 100, -100, 100 };
			sbrtSide.PositionToleranceMm = new List<int>() { 0, 4, 4 };                        // tolerance for the gradient position
			sbrtSide.GradIndexForCoord = 2;                      // index of gradient position to return (zero based index), i.e. the start of the inner wall

			double[] coordBoxLat = new double[2];
			coordBoxLat[0] = GetCoordinates(cooLeft, valHULeft, sbrtSide.GradientHUPerMm, sbrtSide.DistanceInMm, sbrtSide.PositionToleranceMm, sbrtSide.GradIndexForCoord);
			coordBoxLat[1] = GetCoordinates(cooRight, valHURight, sbrtSide.GradientHUPerMm, sbrtSide.DistanceInMm, sbrtSide.PositionToleranceMm, sbrtSide.GradIndexForCoord);
            //coordBoxLat[2] = ((coordBoxRight + coordBoxLeft) / 2);
            if (coordBoxLat[0] !=0 && coordBoxLat[1] != 0)
            {

				MessageBox.Show("found sides");
			}
			return coordBoxLat;

		}

		private double GetSRSLongCoord(PlanSetup plan, VVector dicomPosition, int coordSRSBottom, double coordBoxLeft, double coordBoxRight)
		{
			var image = plan.StructureSet.Image;
			VVector leftFidusStart = dicomPosition;                // only to get the z-coord of the user origo, x and y coord will be reassigned
			VVector rightFidusStart = dicomPosition;               // only to get the z-coord of the user origo, x and y coord will be reassigned
			leftFidusStart.x = coordBoxLeft + 1;                        // start 3 mm in from gradient found in previous step *****************TODO Check for FFS!!!!!!!!!!!!!!!!!!!!!!!
			rightFidusStart.x = coordBoxRight - 1;							// start 1 pixel in right side
			leftFidusStart.y = coordSRSBottom - 91.5;                  // hopefully between fidusles...
			rightFidusStart.y = leftFidusStart.y;
			double stepsFidusLower = 0.5;             //   probably need sub-mm steps to get the fidusle-positions

			VVector leftFidusUpperEnd = leftFidusStart;     // to get the z-coord and x coord, y coord will be reassigned
			VVector leftFidusLowerEnd = leftFidusStart;
			VVector rightFidusUpperEnd = rightFidusStart;
			VVector rightFidusLowerEnd = rightFidusStart;

			int lowerProfileDistance = 40;									// profile length to include all possible fidusles
			int upperProfileDistance = 115;                                 // u gotta love magic numbers and hardcoded values ;)

			leftFidusLowerEnd.y += lowerProfileDistance;                   // distance containing all fidusles determining the Long in 10 cm steps
			leftFidusUpperEnd.y -= upperProfileDistance;                   // endpoint from dimension of box, 
			rightFidusLowerEnd.y += lowerProfileDistance;                  // distance containing all fidusles determining the Long in 10 cm steps
			rightFidusUpperEnd.y -= upperProfileDistance;                  // endpoint from dimension of box, 

			var samplesFidusUpper = (int)Math.Ceiling((leftFidusStart - leftFidusUpperEnd).Length);
			var samplesFidusLower = (int)Math.Ceiling((leftFidusStart - leftFidusLowerEnd).Length / stepsFidusLower);


			// Since the SRS-box walls flexes, the x-coordinate for upper profile may differ from start to end
			// get the max HU in the upper part of the box ( top-most fidusel ) to determine the final x-value for the profile
			var upperLeftFid = leftFidusUpperEnd;
			upperLeftFid.y += 15;							// profile the last 30 mm 
			leftFidusUpperEnd.x = GetMaxHUX(image, upperLeftFid, leftFidusUpperEnd, 6, samplesFidusLower);						// what is the 3 ?????????????????????????????????????
			leftFidusStart.x = GetMaxHUX(image, leftFidusStart, leftFidusLowerEnd, 4, samplesFidusLower);                       // step up 0.2 mm
			leftFidusLowerEnd.x = leftFidusStart.x;																				// Assuming the bottom part of the wall doesn't flex

			var upperRightFid = rightFidusUpperEnd;
			upperRightFid.y += 15;                           // profile the last 30 mm 
			rightFidusUpperEnd.x = GetMaxHUX(image, upperRightFid, rightFidusUpperEnd, -6, samplesFidusLower);
			rightFidusStart.x = GetMaxHUX(image, rightFidusStart, rightFidusLowerEnd, -4, samplesFidusLower);                       // step up 0.2 mm
			rightFidusLowerEnd.x = rightFidusStart.x;


			int numberOfFidusLeft = GetNumberOfFidus(image, leftFidusStart, leftFidusLowerEnd, lowerProfileDistance * 2);
			double fidusLongLeft = GetLongFidus(image, leftFidusStart, leftFidusUpperEnd, upperProfileDistance * 2);
			int numberOfFidusRight = GetNumberOfFidus(image, rightFidusStart, rightFidusLowerEnd, lowerProfileDistance * 2);
			double fidusLongRight = GetLongFidus(image, rightFidusStart, rightFidusUpperEnd, upperProfileDistance * 2);

			// Also need to check the long coordinate above and below (in z-dir) in case its a boundary case where the number of fidusles 
			// steps up. Only neccesary in case of large value for fidusLong or if a discrepancy between the number of fidusles found left and right,
			// or if the long value is not found. +/- 10 mm shift in z-dir is enough to avoid boundary condition *************TODO: check for FFS!!!!!!!!!!!!!
			// TODO: not very DRY...

			double coordSRSLong = 0;
			if (numberOfFidusLeft != numberOfFidusRight || fidusLongLeft > 97 || fidusLongRight > 97 || fidusLongLeft == 0.0 || fidusLongRight ==0)
            {
				int shiftZ = 10;
				leftFidusStart.z += shiftZ;
				leftFidusLowerEnd.z += shiftZ;
				leftFidusUpperEnd.z += shiftZ;
				rightFidusStart.z += shiftZ;
				rightFidusLowerEnd.z += shiftZ;
				rightFidusUpperEnd.z += shiftZ;

				int nOfFidusLeft1 = GetNumberOfFidus(image, leftFidusStart, leftFidusLowerEnd, lowerProfileDistance * 2);
				double fidusLLeft1 = GetLongFidus(image, leftFidusStart, leftFidusUpperEnd, upperProfileDistance * 2);
				int nOfFidusRight1 = GetNumberOfFidus(image, rightFidusStart, rightFidusLowerEnd, lowerProfileDistance * 2);
				double fidusLRight1 = GetLongFidus(image, rightFidusStart, rightFidusUpperEnd, upperProfileDistance * 2);

				leftFidusStart.z -= 2 * shiftZ;
				leftFidusLowerEnd.z -= 2 * shiftZ;
				leftFidusUpperEnd.z -= 2 * shiftZ;
				rightFidusStart.z -= 2 * shiftZ;
				rightFidusLowerEnd.z -= 2 * shiftZ;
				rightFidusUpperEnd.z -= 2 * shiftZ;

				int nOfFidusLeft2 = GetNumberOfFidus(image, leftFidusStart, leftFidusLowerEnd, lowerProfileDistance * 2);
				double fidusLLeft2 = GetLongFidus(image, leftFidusStart, leftFidusUpperEnd, upperProfileDistance * 2);
				int nOfFidusRight2 = GetNumberOfFidus(image, rightFidusStart, rightFidusLowerEnd, lowerProfileDistance * 2);
				double fidusLRight2 = GetLongFidus(image, rightFidusStart, rightFidusUpperEnd, upperProfileDistance * 2);

				// TODO add check to se if 1 and 2 are resonable
				double coordLong1 = (nOfFidusLeft1 + nOfFidusRight1) * 50 + (fidusLLeft1 + fidusLRight1) / 2;
				double coordLong2 = (nOfFidusLeft2 + nOfFidusRight2) * 50 + (fidusLLeft2 + fidusLRight2) / 2;
				//Check if resonable agreement before assigning the final long coordinate as mean value, hard coded values for uncertainty...
				// moved 10 mm in both directions i.e. expected difference in long is 20 mm
                if (Math.Abs(coordLong2 - coordLong1) > 17 && Math.Abs(coordLong2 - coordLong1) < 23)	
				{
					coordSRSLong = (coordLong1 + coordLong2) / 2;
					MessageBox.Show("first " + coordLong1.ToString("0.0") + "\t second " + coordLong2.ToString("0.0"));
				}
                else
                {
					MessageBox.Show("Problem :first " + coordLong1.ToString("0.0") + "\t second " + coordLong2.ToString("0.0"));
				}

			}
            else
            {
				coordSRSLong = (numberOfFidusLeft + numberOfFidusRight) * 50 + (fidusLongLeft + fidusLongRight) / 2;
				MessageBox.Show("Left side " + fidusLongLeft.ToString("0.0") + "\t Right side " + fidusLongRight.ToString("0.0"));
			}

			return coordSRSLong;

		}




		/// <summary>
		/// gets the x-value where maximum HU is found when stepping the y-profile in the direction (Dicom) and range given in steps of 0.1 mm
		/// </summary>
		/// <param name="image"></param>
		/// <param name="fidusStart"></param>
		/// <param name="fidusEnd"></param>
		/// <param name="dirLengthInmm"></param>
		/// <param name="samples"></param>
		/// <returns></returns>
		public static double GetMaxHUX(Image image, VVector fidusStart, VVector fidusEnd, double dirLengthInmm, int samples )
		{
			double newMax = 0.0;
			List<double> HUTemp = new List<double>();
			List<double> cooTemp = new List<double>();
			double finalXValue = 0;
			for (int s = 0; s < 10*Math.Abs(dirLengthInmm); s++)
			{
				fidusStart.x += 0.1* dirLengthInmm/ Math.Abs(dirLengthInmm);  // ugly way to get the direction                     
				fidusEnd.x = fidusStart.x;

				var profFidus = image.GetImageProfile(fidusStart, fidusEnd, new double[samples]);

				for (int i = 0; i < samples; i++)
				{
					HUTemp.Add(profFidus[i].Value);
					cooTemp.Add(profFidus[i].Position.y);
				}
				if (HUTemp.Max() > newMax)
				{
					newMax = HUTemp.Max();
					finalXValue = fidusStart.x;
				}
				HUTemp.Clear();
				cooTemp.Clear();
			}
			//MessageBox.Show("max HU: " + newMax.ToString("0.0"));
			return finalXValue;
		}



		public int GetNumberOfFidus(Image image, VVector fidusStart, VVector fidusEnd, int samples)
		{
			List<double> valHU = new List<double>();
			List<double> coord = new List<double>();
			double findGradientResult;

			var profFidus = image.GetImageProfile(fidusStart, fidusEnd, new double[samples]);

				for (int i = 0; i < samples; i++)
				{
					valHU.Add(profFidus[i].Value);
					coord.Add(profFidus[i].Position.y);
				}
			var fid = new PatternGradient();
			fid.DistanceInMm = new List<double>() { 0, 2 };        // distance between gradients
			fid.GradientHUPerMm = new List<int>() { 100, -100 };    // smallest number of fidusles is one?  
			fid.PositionToleranceMm = new List<int>() { 0, 2 };                         // tolerance for the gradient position, parameter to optimize depending probably of resolution of profile
			fid.GradIndexForCoord = 0;                      // index of gradient position to return, in this case used only as a counter of number of fidusles
															
			findGradientResult = GetCoordinates(coord, valHU, fid.GradientHUPerMm, fid.DistanceInMm, fid.PositionToleranceMm, fid.GradIndexForCoord);
			// keep adding gradient pattern until no more fidusles found
            while (findGradientResult != 0.0)
            {
				fid.DistanceInMm.Add(3);
				fid.GradientHUPerMm.Add(100);
				fid.PositionToleranceMm.Add(2);
				fid.DistanceInMm.Add(2);
				fid.GradientHUPerMm.Add(-100);
				fid.PositionToleranceMm.Add(2);
				fid.GradIndexForCoord++;
				findGradientResult = GetCoordinates(coord, valHU, fid.GradientHUPerMm, fid.DistanceInMm, fid.PositionToleranceMm, fid.GradIndexForCoord);
			}
			return fid.GradIndexForCoord;
		}


		public double GetLongFidus(Image image, VVector fidusStart, VVector fidusEnd, int samples)
		{
			List<double> valHU = new List<double>();
			List<double> coord = new List<double>();
			double findFirstFidus;
			double findSecondFidus;

			var profFidus = image.GetImageProfile(fidusStart, fidusEnd, new double[samples]);

			for (int i = 0; i < samples; i++)
			{
				valHU.Add(profFidus[i].Value);
				coord.Add(profFidus[i].Position.y);
			}
			var fid = new PatternGradient();
			fid.DistanceInMm = new List<double>() { 0, 2 , 49, 2 , 99, 2 };//};        // distance between gradients
			fid.GradientHUPerMm = new List<int>() { 100, -100 , 70, -70 , 100, -100 };//};    // diagonal fidusle have flacker gradient
			fid.PositionToleranceMm = new List<int>() { 2, 3, 105, 4, 105, 3};                        // tolerance for the gradient position, in this case the maximum distance is approx 105 mm
			fid.GradIndexForCoord = 0;                      // index of gradient position to return (zero based index)

			// Finding position of the gradient start is not enough since the long fidusle is diagonal and also changes width depending of the resolution of the image in z-dir, 
			// have to take the mean position before and after.
			double findFirstFidusStart = GetCoordinates(coord, valHU, fid.GradientHUPerMm, fid.DistanceInMm, fid.PositionToleranceMm, fid.GradIndexForCoord);
			fid.GradIndexForCoord = 1;
			double findFirstFidusEnd = GetCoordinates(coord, valHU, fid.GradientHUPerMm, fid.DistanceInMm, fid.PositionToleranceMm, fid.GradIndexForCoord);
			findFirstFidus = (findFirstFidusStart + findFirstFidusEnd) / 2;
			//Find position of second fidus (diagonal)
			fid.GradIndexForCoord = 2;
			double findSecondFidusStart = GetCoordinates(coord, valHU, fid.GradientHUPerMm, fid.DistanceInMm, fid.PositionToleranceMm, fid.GradIndexForCoord);
			fid.GradIndexForCoord = 3;
			double findSecondFidusEnd = GetCoordinates(coord, valHU, fid.GradientHUPerMm, fid.DistanceInMm, fid.PositionToleranceMm, fid.GradIndexForCoord);
			findSecondFidus = (findSecondFidusStart + findSecondFidusEnd) / 2;

			return Math.Abs(findSecondFidus-findFirstFidus);
		}


		/// <summary>
		/// getCoordinates gives the Dicom-coordinates of a gradient 
		/// </summary>
		/// <param name="coord"> 1D coordinates of a profile</param>
		/// <param name="valueHU"> HU-valus of the profile</param>
		/// <param name="hUPerMm"> Gradient to search for in HU/mm with sign indicating direction</param>
		/// <param name="distMm"> Distance in mm to the next gradient</param>
		/// <param name="posTolMm"> Tolerance of position of found gradient in mm</param>
		/// <returns></returns>
		public double GetCoordinates(List<double> coord, List<double> valueHU, List<int> hUPerMm, List<double> distMm, List<int> posTolMm, int indexToReturn)
		{
			string debug = "";
			double[] grad = new double[coord.Count - 1];
			double[] pos = new double[coord.Count - 1];
			int index = 0;

			double gradientStart;
			double gradientEnd;
			double gradientMiddle;
			// resample profile to gradient with position inbetween profile points ( number of samples decreases with one)
			for (int i = 0; i < coord.Count - 2; i++)
			{
				pos[i] = (coord[i] + coord[i + 1]) / 2;
				grad[i] = (valueHU[i + 1] - valueHU[i]) / Math.Abs(coord[i + 1] - coord[i]);
			}



			List<double> gradPosition = new List<double>();
			int indexToReturnToInCaseOfFail = 0;

			for (int i = 0; i < pos.Count(); i++)
			{
				if (index == hUPerMm.Count())                        //break if last condition passed 
				{
					break;
				}
				// if gradient larger than given gradient and in the same direction
				if (Math.Abs((valueHU[i + 1] - valueHU[i]) / Math.Abs(coord[i + 1] - coord[i])) > (Math.Abs(hUPerMm[index])) && SameSign(grad[i], hUPerMm[index]))
				{
					gradientStart = pos[i];
					gradientEnd = pos[i];
					/*/Might not be the largest gradient in the vicinity even if conditions met, step to next and check. If steeper, take that position instead 
					while (Math.Abs((valueHU[i + 2] - valueHU[i + 1])) > Math.Abs((valueHU[i + 1] - valueHU[i])) && SameSign((valueHU[i + 2] - valueHU[i + 1]), (valueHU[i + 1] - valueHU[i])) && i < coord.Count - 1)
					{
						i++;
					}*/
					//Keep stepping up while gradient larger than given huPerMm
					while (Math.Abs((valueHU[i + 1] - valueHU[i]) / Math.Abs(coord[i + 1] - coord[i])) > (Math.Abs(hUPerMm[index])) && SameSign(grad[i], hUPerMm[index]) && i < coord.Count - 2)
					{
						i++;
						gradientEnd = pos[i];
                        if (index == 0)
                        {
							indexToReturnToInCaseOfFail = i + 1; // if the search fails, i.e. can not find next gradient within the distance given, return to position directly after first gradient ends
                        }
					}
					gradientMiddle = (gradientStart + gradientEnd) / 2;
					if (index == 0)     // if this is the first gradient (i.e. index == 0), cannot yet compare the distance between the gradients, step up index and continue
					{
						//gradPosition.Add(pos[i]);
						gradPosition.Add(gradientMiddle);

						index++;
					}
					// if gradient found before expected position (outside tolerance), keep looking
					else if (Math.Abs(gradientMiddle - gradPosition[index - 1]) < distMm[index] - posTolMm[index] && i < pos.Count()-2)
					{
						i++;
						//MessageBox.Show(Math.Abs(gradientMiddle - gradPosition[index - 1]).ToString("0.0"));
					}
					// if next gradient not found within tolerance distance, means that the first gradient is probably wrong, reset index
					else if ((Math.Abs(gradientMiddle - gradPosition[index - 1]) > (Math.Abs(distMm[index]) + posTolMm[index])))
					{
						debug += "Fail " + (Math.Abs(gradientMiddle - gradPosition[index - 1])).ToString("0.0") + "\t" + (distMm[index] + posTolMm[index]).ToString("0.0") + "\n";
						gradPosition.Clear();
						index = 0;
						i = indexToReturnToInCaseOfFail;
					}
					//  compare the distance between the gradients to the criteria given, step up index and continue if within tolerance
					else if ((Math.Abs(gradientMiddle - gradPosition[index - 1]) > (distMm[index] - posTolMm[index])) && (Math.Abs(gradientMiddle - gradPosition[index - 1]) < (distMm[index] + posTolMm[index]))) 
					{
						//debug += pos[i].ToString("0.0") + "\t" + (gradPosition[index] - gradPosition[index - 1]).ToString("0.0") + "\t" + grad[i].ToString("0.0") + "\n";
						gradPosition.Add(gradientMiddle);
						index++;
                        if (index == 1)
                        {
							indexToReturnToInCaseOfFail = i;
                        }
					}
					else
					{   // if not the first gradient and the distance betwen the gradients are not met within the tolerance; reset index and positions and continue search
                        // reset search from second gradient position to avoid missing the actual gradient.
                        if (gradPosition.Count > 1 && indexToReturnToInCaseOfFail > 0)
                        {
							i = indexToReturnToInCaseOfFail;
						}
						gradPosition.Clear();
						index = 0;
					}
				}
			}
			if (index == hUPerMm.Count())
			{

                //debug ***********************
               if (distMm[1] == 4.4)
                {
					MessageBox.Show("Dessa grad bottom: " + Math.Abs(gradPosition[index-3] - gradPosition[index - 4]).ToString("0.0") + "\t" + Math.Abs(gradPosition[index - 2] - gradPosition[index - 3]).ToString("0.0") +
						"\t" +Math.Abs(gradPosition[index - 1] - gradPosition[index - 2]).ToString("0.0"));

				}

				return gradPosition[indexToReturn];
				
			}
			else
			{
				//debug ***********************
				if (distMm[1] == 4.4)
				{
					MessageBox.Show(debug);
				}
				return 0.0;
				
			}
		} // end method GetCoordinates



	}
}
