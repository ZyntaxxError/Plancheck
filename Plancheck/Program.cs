
///////////////////////////////////////////////////////////////////////////////
//  SRTcheck.cs
//
//  ESAPI v15.5 Script for simple plan parameter checks
//  
////////////////////////////////////////////////////////////////////////////////

using System;
using System.Windows;
using System.Text;
using System.Windows.Forms;
using System.Linq;
using VMS.TPS.Common.Model.API;
using VMS.TPS.Common.Model.Types;


// TODO:  Kolla också statiska fält och energi (FFF)
// Couchcheck, VVector[][] GetContoursOnImagePlane eller nåt som säger om rätt bord valt
// DoseDynamic for MLC -> IMRT,  case arc and static MLC...
// TODO: separate the checks for srt
// Naming conventions; structure set, image, reference point, etc
// Maybe just remarks and thats it...? otherwise; no comments!
// VMAT rules;  max 18 cm i x
// Arc rules ;   Om bara en halv- eller trekvartsvarv, bör gå från under pat till över pat.   Om flera arcar, CW och CCW 
// Går ej att kolla om Course diagnos attached... verkar ej heller gå att se vilken algoritm som används för portal dose, bummer...

/// <summary>
/// 
/// </summary>
namespace VMS.TPS         //      Change to VMS.TPS
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
				string courseIntent = context.Course.Intent;
				string courseId = context.Course.Id;


				string message =
				CheckCourseIntent(courseIntent, plan) + "\n" +
				CheckClinProt(plan) +
				CheckCouchStructure(plan.StructureSet) +
				CheckSetupField(plan) + "\n" +
				CheckFieldRules(plan);
				

				MessageBox.Show(message);

				var planIso = plan.Beams.First().IsocenterPosition; // mm från user origo!!!!!!!!! 
				var image = plan.StructureSet.Image;
				var imageUO = image.UserOrigin;             // mm från origo satt från CT!!!!!!
				var imageCTO = image.Origin;                // Ursprungligt origo satt från CT!!!!!! mm från övre vänstra hörnet?

				double imageSizeX = image.XRes * image.XSize;
				//double imageSizeY = image.YRes*image.YSize;
				//double imageSizeZ = image.ZRes*image.ZSize;
				int protocolFraktionIndex = plan.ProtocolID.IndexOf('#');

				MessageBox.Show("Plan iso: " + planIso.x.ToString("0.00") + "\t" + planIso.y.ToString("0.00") + "\n" +
				"User Origo: " + imageUO.x.ToString("0.00") + "\t" + imageUO.y.ToString("0.00") + "\n" +
				"CT origo: " + imageCTO.x.ToString("0.00") + "\t" + imageCTO.y.ToString("0.00") + "\n" +
				"image size x mm :" + imageSizeX + "\n");

				//var imageRes = new double[] {image.XRes,image.YRes,image.ZRes};		// voxel size in mm
				//var imageVoxSize = new int[] {image.XSize, Image.YSize, Image.ZSize};		// image size in voxels










			}
		}


		/*public ( ImageProfile , ImageProfile , ImageProfile ) getImageProfilesThroughIsocenter ( PlanSetup plan )
		{
			var image = plan.StructureSet.Image;
			var dirVecs = new VVector [] {
			plan.StructureSet.Image.XDirection ,
			plan.StructureSet.Image.YDirection ,
			plan.StructureSet.Image.ZDirection
			};
			var steps = new double [] {
			plan.StructureSet.Image.XRes,
			plan.StructureSet.Image.YRes,
			plan.StructureSet.Image.ZRes
			};
			var planIso = plan.Beams.First().IsocenterPosition;
			var tmpRes = new ImageProfile [3];
			//
			// Throws if plan does not have ’BODY ’
			//
			var body = plan.StructureSet.Structures.Single(st => st.Id =="BODY") ;
			for ( int ind = 0; ind < 3; ind ++)
			{
			(var startPoint , var endPoint ) = Helpers.GetStructureEntryAndExit(body ,dirVecs[ind],planIso, steps[ind]) ;
			var samples = (int)Math.Ceiling(( endPoint - startPoint ).Length/steps[ind]) ;
			tmpRes [ind] = image.GetImageProfile( startPoint , endPoint , new double[samples]) ;
			}
		return ( tmpRes [0] , tmpRes [1] , tmpRes [2]) ;
		}*/





		// ********* 	Kontroll om SRT-plan (enbart baserad på fraktionering, aning klent och funkar ej...) *********

		public bool IsPlanSRT(PlanSetup plan)
		{
			return (plan.NumberOfFractions > 2 && plan.DosePerFraction.Dose > 6.0 && plan.TotalDose.Dose >= 45.0);
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


		// ********* 	Kontroll av kliniskt protokoll-ID om sådant kopplat till plan, jämförs med plan-fraktionering *********

		public string CheckClinProt(PlanSetup plan)
		{
			string cResults = "";
			if (plan.ProtocolID.Length != 0)
			{
				int protocolFractionIndex = plan.ProtocolID.IndexOf('#');										// find the index of the symbol indicating nr of fractions
				string protocolFrNrInfo = plan.ProtocolID.Substring(protocolFractionIndex - 2, 2).Trim();       // retrieve the two characters before the #, and remove whitespaces
				if (Int32.TryParse(protocolFrNrInfo, out int fractionsInProtocol))								// try parsing it to int and, if successful, compare to plan fractions
				{
					if (fractionsInProtocol != plan.NumberOfFractions)
					{
						cResults = "** Check the attached clinical protocol! \n \n";
					}
				}
			}
			return cResults;
		}

		// ********* 	Targetvolym, referenspunkt *********

		public string CheckPlanProp(PlanSetup plan, StructureSet sSet)
		{
			string cResults = "";
			/*var target = (from s in sSet.Structures
						  where s.Id == plan.TargetVolumeID
						  select s);*/

			if (string.IsNullOrEmpty(plan.TargetVolumeID))
			{
				cResults = "** No plan target volume selected \n";
			}
			else
			{
				Structure target = sSet.Structures.Where(s => s.Id == plan.TargetVolumeID).Where(s => s.DicomType == "PTV").FirstOrDefault();
				if (target == null)
				{
					cResults = "** Plan target volume should be of type PTV \n";
				}
				else
				{
					cResults = "Plan target volume: " + target.Id;
				}
			}
			return cResults;
		}

			   		 	  	  


		// ********* 	Kontroll av att bordsstruktur existerar, inte är tom och har korrekt HU 	********* 
		// begränsningar: kollar ej att rätt bordstyp används och ej heller positionering i förhållande till body

		public string CheckCouchStructure(StructureSet SSet)
		{
			string cResult = "";
			bool couchExt = false;
			bool couchInt = false;
			foreach (Structure s in SSet.Structures)
			{
				if (s.Id.Contains("CouchSurf") && !s.IsEmpty)
				{
					s.GetAssignedHU(out double couchExtHU);
					if (Math.Round(couchExtHU) == -300)
					{
						couchExt = true;
					}
				}
				if (s.Id.Contains("CouchInt") && !s.IsEmpty)
				{
					s.GetAssignedHU(out double couchIntHU);
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



		// ********* 	Kontroll av Setup-fält, namngivning	********* 

		public string CheckSetupField(PlanSetup plan)
		{
			string cResults = "";
			int countSetupfields = 0;
			foreach (var beam in plan.Beams)
			{
				if (beam.IsSetupField)
				{
					countSetupfields++;
					if (beam.Id.ToUpper().Substring(0, 2).Equals(plan.Id.ToUpper().Substring(0, 2)))
					{
						if (!beam.Id.ToUpper().Contains("CBCT"))
						{
							cResults = cResults + CheckPlanarSetupFields(beam) + "\n";
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



		public string CheckPlanarSetupFields(Beam beam)
		{
			string cResults = "";
			string trimmedID = beam.Id.Substring(2).Trim();         // start iteration at index 2 (PX index 0 and 1)
			int gantryAngleInBeamID = 1000;
			for (int i = 1; i <= trimmedID.Length; i++)
			{
				if (Int32.TryParse(trimmedID.Substring(0, i), out int test))
				{
					gantryAngleInBeamID = test;
				}
			}
			if (gantryAngleInBeamID != 1000)
			{
				if (gantryAngleInBeamID == Math.Round(beam.ControlPoints.First().GantryAngle))
				{
					cResults = "OK";
				}
				else
				{
					cResults = "Check setup field ID!";
				}
			}
			return cResults;
		}

			   		 


		// ********* 	Kontroll av diverse fältregler och "best practices"	********* 


		public string CheckFieldRules(PlanSetup plan)
		{
			string remarksSRS = "";
			string remarksDoseRateFFF = "";
			string remarksArcDynFFF = "";
			string remarksArcDynCollAngle = "";
			string remarksArcStartStop = "";
			string remarksArcDirections = "";

			int countFields = plan.Beams.Count();
			//	int countWedge = plan.Beams.Wedges.Count();
			int countSetupFields = plan.Beams.Where(b => b.IsSetupField).Count();
			int countTreatFields = countFields - countSetupFields;
			int countRemarksSRS = 0;                    // All fields should be SRS if SBRT- or SRS-plan
			int countRemarksDoseRateFFF = 0;            // Dose rate should be maximum for FFF
			int countRemarksArcDynFFF = 0;              // The energy should be FFF if dynamic arc used
			int countRemarksArcDynCollAngle = 0;        // The collimator angle should be between +/-5 deg if dynamic arc used
			int countArcCW = 0;
			int countArcCCW = 0;                        // the absolute difference between CW and CCW should be less than two...
			foreach (var beam in plan.Beams)
			{
				if (!beam.IsSetupField)
				{
					if (!beam.Technique.Id.Contains("SRS") && IsPlanSRT(plan))
					{
						if (countRemarksSRS < 1) { remarksSRS = "** Change technique to SRS-" + beam.Technique.Id + "! \n"; };
						countRemarksSRS++;
					}
					if (beam.EnergyModeDisplayName.Contains("FFF"))
					{
						remarksDoseRateFFF += CheckDoseRateFFF(beam, ref countRemarksDoseRateFFF);
					}
					if (beam.MLCPlanType == MLCPlanType.ArcDynamic && IsPlanSRT(plan))
					{
						remarksArcDynFFF += CheckArcDynFFF(beam, ref countRemarksArcDynFFF);
						remarksArcDynCollAngle += CheckArcDynCollAngle(beam, ref countRemarksArcDynCollAngle);
					}
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
						remarksArcStartStop += CheckArcStartStop(beam);
					}
				}

			}
			if (Math.Abs(countArcCCW - countArcCW) > 1)
			{
				remarksArcDirections += "** Check the arc directions! \t";
			}
			string cResults = remarksSRS + remarksDoseRateFFF + remarksArcDynFFF + remarksArcDynCollAngle + remarksArcStartStop + remarksArcDirections;
			return cResults;
		}


		// ********* 	Kontroll av dosrat vid FFF	*********

		public string CheckDoseRateFFF(Beam beam, ref int countRemarksDoseRateFFF)
		{
			string cResults = "";
			if (beam.DoseRate < 1400)
			{
				if (countRemarksDoseRateFFF < 1)
				{
					cResults = "** Change dose rate to maximum! Field: " + beam.Id;
				}
				else
				{
					cResults += ", " + beam.Id;
				}
				countRemarksDoseRateFFF++;
			}
			return cResults;
		}


		// ********* 	Kontroll av energi vid Dynamic Arc	*********

		public string CheckArcDynFFF(Beam beam, ref int countRemarksArcDynFFF)
		{
			string cResults = "";
			if (countRemarksArcDynFFF < 1 && !beam.EnergyModeDisplayName.Contains("FFF"))
			{
				cResults = "** Change energy to FFF! \n";
				countRemarksArcDynFFF++;
			}
			return cResults;
		}


		// ********* 	Kontroll av kollimatorvinkel vid Dynamic Arc	*********

		public string CheckArcDynCollAngle(Beam beam, ref int countRemarksArcDynCollAngle)
		{
			string cResults = "";
			if (countRemarksArcDynCollAngle < 1 && beam.ControlPoints.First().CollimatorAngle > 5.0 && beam.ControlPoints.First().CollimatorAngle < 355.0)
			{
				cResults = "** Collimator angle for DynArc should be between +/- 5 deg \n";
				countRemarksArcDynCollAngle++;
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
