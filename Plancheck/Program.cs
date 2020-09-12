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
				string courseIntent = context.Course.Intent;
				string courseId = context.Course.Id;

				string message =
				CheckCourseIntent(courseIntent, plan) + "\n" +
				CheckClinProt(plan) +
				CheckPlanProp(plan, plan.StructureSet) +
				CheckCouchStructure(plan.StructureSet) +
				CheckSetupField(plan) + "\n" +
				"Treatment fields \n" +
				"ID \t Energy \t Tech. \t Drate \t MLC \n" +
				CheckFieldRules(plan);

				MessageBox.Show(message);



				// testing ground





				var planIso = plan.Beams.First().IsocenterPosition; // mm from Dicom-origo
				var image = plan.StructureSet.Image;
				var imageUserOrigo = image.UserOrigin;             // mm från origo satt från CT vilket är dicom-origo!The user origin in DICOM coordinates in millimeter. 
				var imageCTO = image.Origin;                // Ursprungligt origo satt från CT, mm från CT-origo TILL övre vänstra hörnet i första bilden!
															//The origin of the image. In other words, the DICOM coordinates of the center point of the upper-left hand corner voxel of the first image plane

				double imageSizeX = image.XRes * image.XSize;
				double imageSizeY = image.YRes * image.YSize;
				double dist = VVector.Distance(imageCTO, imageUserOrigo);  // 3D distance from one koord to another, double 
				var userIsoCoord = image.DicomToUser(planIso, plan);  // 

				// VVector @ image center in iso plane in dicomcoordinates
				VVector isoPlaneImageCenter = planIso;




				double xLeftUpperCorner = image.Origin.x - image.XRes / 2;  // Dicomcoord in upper left corner ( NOT middle of voxel in upper left corner)
				double yLeftUpperCorner = image.Origin.y - image.YRes / 2;  // Dicomcoord in upper left corner ( NOT middle of voxel in upper left corner)


				// instead of absolute image center => coord of first and last image voxel in x and y, XSize in voxels and XRes in mm/voxel
				// this seems like a smart idea if I can understand it...
				double xVoxStart = image.Origin.x;
				double xVoxEnd = image.Origin.x + image.XRes * image.XSize - image.XRes;
				double yVoxStart = image.Origin.y;
				double yVoxEnd = image.Origin.y + image.YRes * image.YSize - image.YRes;



				//**********  Image profile, need to define startpoint(VVector), endpoint(VVector) and number of samples (double[]) in that distance *********** 


				//In case of SBRT we want to get a profile in User origo plane, middle of the image in x-direction and starting from bottom of the image in y.


				VVector bottomProfileStart = image.UserOrigin;              // only to get the z-coord of the user origo, x and y coord will be reassigned
				bottomProfileStart.x = xLeftUpperCorner + imageSizeX / 2;           // center of the image in x-direction
				bottomProfileStart.y = yLeftUpperCorner + imageSizeY - image.YRes;  // start 1 pixel in from bottom...
				double steps = plan.StructureSet.Image.YRes;                //   (mm/voxel) to make the steps 1 pixel wide, can skip this if 1 mm steps is wanted

				VVector bottomProfileEnd = bottomProfileStart;
				bottomProfileEnd.y -= 200 * steps;                  // endpoint 200 steps in -y direction

				var samplesY = (int)Math.Ceiling((bottomProfileStart - bottomProfileEnd).Length / steps);
				var profY = image.GetImageProfile(bottomProfileStart, bottomProfileEnd, new double[samplesY]);

				/*
							MessageBox.Show("Plan iso: " + planIso.x.ToString("0.0") + "\t" + planIso.y.ToString("0.0") + "\n" +
							"User Origo: " + imageUserOrigo.x.ToString("0.0")  + "\t" + imageUserOrigo.y.ToString("0.0")  + "\n" +
							"CT origo: " + imageCTO.x.ToString("0.0") + "\t" + imageCTO.y.ToString("0.0")  +  "\n" +
							"image size x mm :" + imageSizeX + "\t" + "y: " + imageSizeY + "\n" +
							"DicToUserIso: " + userIsoCoord.x.ToString("0.00") +  "\t" + userIsoCoord.y.ToString("0.00") +  "\n" +
							userIsoCoord.y.ToString("0.00") +  "\n" +
							//isoPlaneImageCenter.y
							//profY[1].Position.x +  "\n" +
							//profY[1].Value.ToString("0.00") +  "\n" +					// seems to give value directly in HU
							//image.VoxelToDisplayValue(Convert.ToInt32(profY[1].Value)) +  "\n" +
							"Number of samples\t" + samplesY + "\n" +
							//"PlaneImageCenter\t" + isoPlaneImageCenter.y + "\n" +
							image.XDirection.y);
				*/

				//var imageRes = new double[] {image.XRes,image.YRes,image.ZRes};		// voxel size in mm
				//var imageVoxSize = new int[] {image.XSize, Image.YSize, Image.ZSize};		// image size in voxels


				//***********  Gradient patter describing expected profile in HU of the Lax-box bottom **********

				PatternGradient lax = new PatternGradient();
				lax.DistanceInMm = new List<double>() { 0, 4.4, 12.3, 4.4 };        // distance between gradients, mean values from profiling 10 pat 
				lax.GradientHUPerMm = new List<int>() { 100, -100, 100, -100 };
				lax.PositionToleranceMm = 2;                        // tolerance for the gradient position
				lax.gradIndexForCoord = 2;                      // index of gradient position to return (zero based index)



				// Imageprofile gets a VVector back, take the coordinates and respective HU and put them in two Lists of double, might be better ways of doing this...
				List<double> valHU = new List<double>();
				List<double> coo = new List<double>();
				string debug1 = "";
				for (int i = 0; i < samplesY; i++)
				{
					valHU.Add(profY[i].Value);
					coo.Add(profY[i].Position.y);
					debug1 += coo[i].ToString("0.0") + "\t" + valHU[i].ToString("0.0") + "\n";
				}
				//MessageBox.Show(debug1);

				// Get the coordinate (dicom) that represents inner bottom of laxbox (-2 in box-coordinates)
				double coordBoxBottom = getCoordinates(coo, valHU, lax.GradientHUPerMm, lax.DistanceInMm, lax.PositionToleranceMm, lax.gradIndexForCoord);
				// in Boxcoordinates this is equal to -2 in ant-post, can then check the coordinates for user origo in y which should be 95 (in SRS coord) by adding -97 to found coordinate

				if (coordBoxBottom != 0)
				{

					// ************************************ get profiles in x direction, left and right side and detemine center of box ********************



					VVector leftProfileStart = image.UserOrigin;                // only to get the z-coord of the user origo, x and y coord will be reassigned
					VVector rightProfileStart = image.UserOrigin;               // only to get the z-coord of the user origo, x and y coord will be reassigned
					leftProfileStart.x = xLeftUpperCorner + image.XRes;         // start 1 pixel in left side
					rightProfileStart.x = xLeftUpperCorner + image.XSize * image.XRes - image.XRes;         // start 1 pixel in right side
					leftProfileStart.y = coordBoxBottom - 93.5;                 // hopefully between fidusles...
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

					PatternGradient laxSide = new PatternGradient();
					laxSide.DistanceInMm = new List<double>() { 0, 2, 13 };     // distance between gradients, mean values from profiling 10 pat 
					laxSide.GradientHUPerMm = new List<int>() { 100, -100, 100 };
					laxSide.PositionToleranceMm = 2;                        // tolerance for the gradient position
					laxSide.gradIndexForCoord = 2;                      // index of gradient position to return (zero based index)


					double coordBoxLeft = getCoordinates(cooLeft, valHULeft, laxSide.GradientHUPerMm, laxSide.DistanceInMm, laxSide.PositionToleranceMm, laxSide.gradIndexForCoord);
					double coordBoxRight = getCoordinates(cooRight, valHURight, laxSide.GradientHUPerMm, laxSide.DistanceInMm, laxSide.PositionToleranceMm, laxSide.gradIndexForCoord);

					//MessageBox.Show(coordBoxLeft.ToString("0.0") + "\t" + coordBoxRight.ToString("0.0") + "\n" + ((coordBoxRight+coordBoxLeft)/2).ToString("0.0") );


					string UserOrigoCheck = "";
					if (Math.Abs(image.UserOrigin.y - (coordBoxBottom - 96)) < 3 && Math.Abs(image.UserOrigin.x - ((coordBoxRight + coordBoxLeft) / 2)) < 3)
					{
						UserOrigoCheck = "User origo position in SRS frame seems OK in Lat and Vrt." + "\n\n" + "Estimated position in SRS frame coordinates from image profiles: \n\n" +
						" Lat: " + (-(image.UserOrigin.x - ((coordBoxRight + coordBoxLeft) / 2) - 300)).ToString("0.0") + "\t Vrt: " + Math.Abs(image.UserOrigin.y - (coordBoxBottom - 1)).ToString("0.0") + "\t (+/- 2mm)";
					}
					else if (coordBoxRight == 0 || coordBoxLeft == 0)
					{
						UserOrigoCheck = "Cannot find the SRS-frame, no automatic check of User origo possible.";
					}
					else
					{
						UserOrigoCheck = "Check position of user origo";
					}
					MessageBox.Show(UserOrigoCheck);

				}





				// TODO check if isocenter in same plane as user origo, not neccesary though as there can be multiple isocenters (muliple plans)


				//*********************************    "normal" checks *********************************




			}
		}



		// ********* 	Kontroll om SRT-plan (enbart baserad på fraktionering, aning klent men bör fungera) *********

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
			int fractionsInProtocol = 0;
			if (plan.ProtocolID.Length != 0)
			{
				int protocolFractionIndex = plan.ProtocolID.IndexOf('#');                   // find the index of the symbol indicating nr of fractions
				string protocolFrNrInfo = plan.ProtocolID.Substring(protocolFractionIndex - 2, 2).Trim();       // retrieve the two characters before the #, and remove whitespaces
				if (Int32.TryParse(protocolFrNrInfo, out fractionsInProtocol))                  // try parsing it to int and, if successful, compare to plan fractions
				{
					if (fractionsInProtocol != plan.NumberOfFractions)
					{
						cResults = "** Check the attached clinical protocol! \n \n";
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
					cResults = "* Plan target volume should be of type PTV \n";
				}
				else
				{
					cResults = "Plan target volume: " + target.Id;
				}
			}
			return cResults + "\n";
		}





		// ********* 	Kontroll av att bordsstruktur existerar, är av rätt typ, inte är tom och har korrekt HU 	********* 
		// begränsningar: kollar ej positionering i förhållande till body, rätt bordstyp kollas enbart på namn
		// Exact IGRT Couch, thin, medium och thick kan vara korrekt beroende på lokalisation...

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
					if (beam.Id.ToUpper().Substring(0, 2).Equals(plan.Id.ToUpper().Substring(0, 2)))
					{
						if (beam.Id.ToUpper().Contains("CBCT"))
						{
							cResults = cResults + "OK" + "\n";
						}
						else
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

		// ********* 	Kontroll av planara Setup-fält (ej namngivna cbct), namngivning efter gantryvinkel	********* 

		public string CheckPlanarSetupFields(Beam beam)
		{
			string cResults = "";
			string trimmedID = beam.Id.Substring(2).Trim();         // start iteration at index 2 (PX index 0 and 1)
			int gantryAngleInBeamID = 1000;
			int test = 0;
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
					cResults += "OK";
				}
				else
				{
					cResults += "Check name!";
				}
			}
			else
			{
				cResults += "Check name!";
			}
			return cResults;
		}



		// ********* 	Kontroll av diverse fältregler och "best practices"	********* 


		public string CheckFieldRules(PlanSetup plan)
		{
			string cResults = "";
			string remarks = "";
			int countFields = plan.Beams.Count();
			int countSetupFields = plan.Beams.Where(b => b.IsSetupField).Count();
			int countTreatFields = countFields - countSetupFields;
			int countSRSRemarks = 0;                    // All fields should be SRS if SBRT- or SRS-plan
			int countDoseRateFFFRemarks = 0;            // Dose rate should be maximum for FFF
			int countArcDynFFFRemarks = 0;              // The energy should be FFF if dynamic arc used
			int countArcDynCollAngleRemarks = 0;        // The collimator angle should be between +/-5 deg if dynamic arc used
			int countArcCW = 0;
			int countArcCCW = 0;                        // the absolute difference between CW and CCW should be less than two...
			foreach (var beam in plan.Beams)
			{
				if (!beam.IsSetupField)
				{
					cResults = cResults + beam.Id + "\t" + beam.EnergyModeDisplayName + "\t" + beam.Technique.Id + "\t" + beam.DoseRate + "\t" + beam.MLCPlanType + "\n";
					if (!beam.Technique.Id.Contains("SRS"))
					{
						if (countSRSRemarks < 1) { remarks = remarks + "** Change technique to SRS-" + beam.Technique.Id + "! \n"; };
						countSRSRemarks++;
					}
					if (beam.EnergyModeDisplayName.Contains("FFF"))
					{
						remarks += CheckDoseRateFFF(beam, ref countDoseRateFFFRemarks);
					}
					if (beam.MLCPlanType == MLCPlanType.ArcDynamic)
					{
						remarks += CheckArcDynFFF(beam, ref countArcDynFFFRemarks);
						remarks += CheckArcDynCollAngle(beam, ref countArcDynCollAngleRemarks);
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
						remarks += CheckArcStartStop(beam);
					}
				}
				if (Math.Abs(countArcCCW - countArcCW) > 1)
				{
					remarks += "** Check the arc directions! \t";
				}
			}
			return cResults + "\n" + remarks;
		}


		// ********* 	Kontroll av dosrat vid FFF	*********

		public string CheckDoseRateFFF(Beam beam, ref int countDoseRateFFFRemarks)
		{
			string cResults = "";
			if (beam.DoseRate < 1400)
			{
				if (countDoseRateFFFRemarks < 1)
				{
					cResults = "** Change dose rate to maximum! \n";
					countDoseRateFFFRemarks++;
				}
			}
			return cResults;
		}


		// ********* 	Kontroll av energi vid Dynamic Arc	*********

		public string CheckArcDynFFF(Beam beam, ref int countArcDynFFFRemarks)
		{
			string cResults = "";
			if (countArcDynFFFRemarks < 1 && !beam.EnergyModeDisplayName.Contains("FFF"))
			{
				cResults = "** Change energy to FFF! \n";
				countArcDynFFFRemarks++;
			}
			return cResults;
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

		bool sameSign(double num1, double num2)
		{
			return num1 >= 0 && num2 >= 0 || num1 < 0 && num2 < 0;
		}

		public class PatternGradient
		{
			public List<double> DistanceInMm { get; set; }
			public List<int> GradientHUPerMm { get; set; }
			public int PositionToleranceMm { get; set; }
			public int gradIndexForCoord { get; set; }
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
		double getCoordinates(List<double> coord, List<double> valueHU, List<int> hUPerMm, List<double> distMm, int posTolMm, int indexToReturn)
		{
			string debug = "";
			double[] grad = new double[coord.Count - 1];
			double[] pos = new double[coord.Count - 1];
			int index = 0;

			for (int i = 0; i < coord.Count - 1; i++)
			{
				pos[i] = (coord[i] + coord[i + 1]) / 2;
				grad[i] = (valueHU[i + 1] - valueHU[i]) / Math.Abs(coord[i + 1] - coord[i]);
			}


			List<double> gradPosition = new List<double>();

			for (int i = 0; i < coord.Count - 1; i++)
			{
				pos[i] = (coord[i] + coord[i + 1]) / 2;
				grad[i] = (valueHU[i + 1] - valueHU[i]) / Math.Abs(coord[i + 1] - coord[i]);
				if (index == hUPerMm.Count())                        //break if last condition passed 
				{
					break;
				}
				// if gradient larger than given gradient and in the same direction
				if (Math.Abs((valueHU[i + 1] - valueHU[i]) / Math.Abs(coord[i + 1] - coord[i])) > (Math.Abs(hUPerMm[index])) && sameSign(grad[i], hUPerMm[index]))
				{
					//Might not be the largest gradient in the visinity even if conditions met, step to next and check. If steeper, take that position instead 
					while (Math.Abs((valueHU[i + 2] - valueHU[i + 1])) > Math.Abs((valueHU[i + 1] - valueHU[i])) && sameSign((valueHU[i + 2] - valueHU[i + 1]), (valueHU[i + 1] - valueHU[i])))
					{
						i++;
					}
					gradPosition.Add(pos[i]);
					if (index == 0)     // if this is the first gradient (i.e. index == 0), cannot yet compare the distance between the gradients, step up index and continue
					{
						debug += pos[i].ToString("0.0") + "\t" + grad[i].ToString("0.0") + "\n";
						index++;
					}
					//  compare the distance between the gradients to the criteria given, step up index and continue
					else if ((Math.Abs(gradPosition[index] - gradPosition[index - 1]) > (distMm[index] - posTolMm)) && (Math.Abs(gradPosition[index] - gradPosition[index - 1]) < (distMm[index] + posTolMm))) // jämför avstånd mellan gradienter mot angett avstånd +/- marginal
					{
						debug += pos[i].ToString("0.0") + "\t" + (gradPosition[index] - gradPosition[index - 1]).ToString("0.0") + "\t" + grad[i].ToString("0.0") + "\n";
						//gradPosition.Add(pos[i]);
						index++;
					}
					else
					{                   // if not first gradient and distance betwen the gradients not met within the tolerance, reset index and positions and continue search
						gradPosition.Clear();
						index = 0;
					}
				}
			}
			if (index == hUPerMm.Count())
			{
				return gradPosition[indexToReturn];
			}
			else
			{
				return 0;
			}
		} // end method GetCoordinates



	}
}
