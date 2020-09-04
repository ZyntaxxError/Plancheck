////////////////////////////////////////////////////////////////////////////////
//  SRTcheck.cs
//
//  ESAPI v15.5 Script for simple plan parameter checks
//  
////////////////////////////////////////////////////////////////////////////////

using System;
using System.Windows;
using System.Text;
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





				var planIso = plan.Beams.First().IsocenterPosition; // mm från user origo!!!!!!!!!? 
				var image = plan.StructureSet.Image;
				var imageUserOrigo = image.UserOrigin;             // mm från origo satt från CT vilket är dicom-origo!!!!!!The user origin in DICOM coordinates in millimeter. 
				var imageCTO = image.Origin;                // Ursprungligt origo satt från CT!!!!!! mm från CT-origo TILL övre vänstra hörnet i första bilden!?
				//The origin of the image. In other words, the DICOM coordinates of the center point of the upper-left hand corner voxel of the first image plane

				double imageSizeX = image.XRes * image.XSize;
				double imageSizeY = image.YRes * image.YSize;
				double dist = VVector.Distance(imageCTO,imageUserOrigo);  // 3D distance from one koord to another, double 
				var userIsoCoord = image.DicomToUser(planIso, plan);  // 

				// VVector @ image center in iso plane in dicomcoordinates
				VVector isoPlaneImageCenter = planIso;

				double xLeftUpperCorner = image.Origin.x-image.XRes/2;	// Dicomcoord in upper left corner ( NOT middle of voxel in upper left corner)
				double yLeftUpperCorner = image.Origin.y+(imageSizeY-image.YRes)/2;	// Dicomcoord in upper left corner ( NOT middle of voxel in upper left corner)
				isoPlaneImageCenter.x = xLeftUpperCorner;
				isoPlaneImageCenter.y = yLeftUpperCorner;

				// instead of absolute image center => coord of first and last image voxel in x and y, XSize in voxels and XRes in mm/voxel
				double xVoxStart = image.Origin.x;
				double xVoxEnd = image.Origin.x + image.XRes * image.XSize - image.XRes;
				double yVoxStart = image.Origin.y;
				double yVoxEnd = image.Origin.y + image.YRes * image.YSize - image.YRes;
				




			double steps =  plan.StructureSet.Image.XRes;
			var endPoint = imageUserOrigo;
			endPoint.x += 5*steps;
			//var profX = new ImageProfile;

			var samples = (int)Math.Ceiling(( endPoint - imageUserOrigo ).Length/steps) ;
			var profX = image.GetImageProfile( planIso , endPoint , new double [samples]) ;


				MessageBox.Show("Plan iso: " + planIso.x.ToString("0.00") + "\t" + planIso.y.ToString("0.00") + "\n" +
				"User Origo: " + imageUserOrigo.x.ToString("0.00")  + "\n" +
				"CT origo: " + imageCTO.x.ToString("0.00") +  "\n" +
				"image size x mm :" + imageSizeX + "\n" +
				userIsoCoord.x.ToString("0.00") +  "\n" +
				userIsoCoord.y.ToString("0.00") +  "\n" +
				profX[1].Position.x +  "\n" +
				Convert.ToInt32(profX[1].Value) +  "\n" +					// seems to give value directly in HU
				image.VoxelToDisplayValue(Convert.ToInt32(profX[1].Value)) +  "\n" +
				"Number of samples\t" + samples + "\n" +
				"PlaneImageCenter\t" + isoPlaneImageCenter.y + "\n" +
				image.XDirection.x);
				

				//var imageRes = new double[] {image.XRes,image.YRes,image.ZRes};		// voxel size in mm
				//var imageVoxSize = new int[] {image.XSize, Image.YSize, Image.ZSize};		// image size in voxels


				            PatternGradient lax = new PatternGradient();
            lax.DistanceInMm = new List<double>();
            lax.GradientHUPerMm = new List<int>();
            lax.DistanceInMm.Add(0);
            lax.GradientHUPerMm.Add(150);
            lax.DistanceInMm.Add(5);
            lax.GradientHUPerMm.Add(-150);
            lax.DistanceInMm.Add(13);
            lax.GradientHUPerMm.Add(150);
            lax.DistanceInMm.Add(5);
            lax.GradientHUPerMm.Add(-150);
				
			//	getCoordinates(List<double> coord, List<double> valueHU);

/*

		public ImageProfile getImageProfileXThroughIsocenter ( PlanSetup plan )
		{

			double steps =  plan.StructureSet.Image.XRes;


			var endPoint = imageUserOrigo;
			endPoint.x += 5;
			var profX = new ImageProfile[1];




			//var samples = (int)Math.Ceiling(( endPoint - startPoint ).Length/steps) ;

			profX = image.GetImageProfile( planIso , endPoint , new double[samples]) ;

		return profX ;
		}

*/










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















			}
		}















		// ********* 	Kontroll om SRT-plan (enbart baserad på fraktionering, aning klent men bör fungera) *********
		
		public bool IsPlanSRT(PlanSetup plan)
			{
			return (plan.NumberOfFractions > 2 && plan.DosePerFraction.Dose >6.0 && plan.TotalDose.Dose >= 45.0); 
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
		if(plan.ProtocolID.Length != 0)
			{
			int protocolFractionIndex = plan.ProtocolID.IndexOf('#');					// find the index of the symbol indicating nr of fractions
			string protocolFrNrInfo = plan.ProtocolID.Substring(protocolFractionIndex-2,2).Trim();		// retrieve the two characters before the #, and remove whitespaces
			if (Int32.TryParse(protocolFrNrInfo, out fractionsInProtocol))					// try parsing it to int and, if successful, compare to plan fractions
				{
				if(fractionsInProtocol != plan.NumberOfFractions)
					{
						cResults = "** Check the attached clinical protocol! \n \n";
					}
				}
			}
			return cResults;
		}

		// ********* 	Targetvolym; kollar att det är valt och av typen PTV *********

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


		// ********* 	Kontroll av Setup-fält, namngivning	********* 

		public string CheckSetupField(PlanSetup plan)
		{
			string cResults = "Setup-field: \t";
			int countSetupfields = 0;
			foreach (var beam in plan.Beams)
			{
				if (beam.IsSetupField)
				{
					cResults = cResults + beam.Id + "\t \t";
					countSetupfields++;
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
		string trimmedID = beam.Id.Substring(2).Trim();			// start iteration at index 2 (PX index 0 and 1)
		int gantryAngleInBeamID = 1000;
		int test = 0;
			for (int i=1; i<trimmedID.Length; i++)				
			{
				if (Int32.TryParse(trimmedID.Substring(0,i), out test))					// try parsing it to int and
				{
					gantryAngleInBeamID = test;
				}
			}
			if(gantryAngleInBeamID != 1000)
			{
				if(gantryAngleInBeamID == Math.Round(beam.ControlPoints.First().GantryAngle))
				{
				cResults = "OK";
				}
				else
				{
				cResults = "Check name!";
				}
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
			int countSRSRemarks = 0;					// All fields should be SRS if SBRT- or SRS-plan
			int countDoseRateFFFRemarks = 0;			// Dose rate should be maximum for FFF
			int countArcDynFFFRemarks = 0;				// The energy should be FFF if dynamic arc used
			int countArcDynCollAngleRemarks = 0;        // The collimator angle should be between +/-5 deg if dynamic arc used
			int countArcCW = 0;
			int countArcCCW = 0;						// the absolute difference between CW and CCW should be less than two...
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
				if (Math.Abs(countArcCCW-countArcCW)>1)
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
        }
		
		
		
                        void getCoordinates(List<double> coord, List<double> valueHU)
            {


                double[] grad = new double[coord.Count - 1];
                double[] pos = new double[coord.Count - 1];
                int index = 0;

                for (int i = 0; i < coord.Count - 1; i++)
                {
                    pos[i] = (coord[i] + coord[i + 1]) / 2;
                    grad[i] = (valueHU[i + 1] - valueHU[i]) / Math.Abs(coord[i + 1] - coord[i]);
                    System.Console.WriteLine(pos[i] + "\t" + grad[i]);
                }
                                                          
                List<double> gradPosition = new List<double>();
                               
                for (int i = 0; i < coord.Count - 1; i++)
                {
                    pos[i] = (coord[i] + coord[i + 1]) / 2;
                    grad[i] = (valueHU[i + 1] - valueHU[i]) / Math.Abs(coord[i + 1] - coord[i]);
                    if (index > lax.GradientHUPerMm.Count() - 1)                        //break if last condition passed 
                    {
                        break;
                    }

                    if (Math.Abs((valueHU[i + 1] - valueHU[i]) / Math.Abs(coord[i + 1] - coord[i])) > (Math.Abs(lax.GradientHUPerMm[index])) && sameSign(grad[i], lax.GradientHUPerMm[index])) // om gradient > angiven gradient och om går åt samma håll
                    {
                        //Might not be the largest gradient in the visinity 
                        while (Math.Abs((valueHU[i + 2] - valueHU[i + 1])) > Math.Abs((valueHU[i + 1] - valueHU[i])) && sameSign((valueHU[i + 2] - valueHU[i + 1]), (valueHU[i + 1] - valueHU[i])))   // THIS IS PROBABLY WRONG
                        {
                            i++;
                        }

                        System.Console.WriteLine(pos[i] + "\t" + grad[i]);

                        gradPosition.Add(pos[i]);

                        if (index == 0)
                        {
                            index++;
                            // stega upp i objektet till nästa gradient om distance är uppfyllt, tills sista är uppfylld och då break
                        }
                        else if ((Math.Abs(gradPosition[index] - gradPosition[index - 1]) > (lax.DistanceInMm[index] - 3)) && (Math.Abs(gradPosition[index] - gradPosition[index - 1]) < (lax.DistanceInMm[index] + 3))) // jämför avstånd mellan gradienter mot angett avstånd +/- marginal
                        {
                            //gradPosition.Add(pos[i]);
                            index++;
                        }
                        else
                        {
                            gradPosition.Clear();       // om det inte är första gradienten och om inte avståndet stämmer; nollställ o börja om
                            index = 0;
                        }

                        foreach (var item in gradPosition)
                        {
                            System.Console.WriteLine(item);
                        }
                        System.Console.WriteLine("\n");
                    }
                }
                foreach (var item in gradPosition)
                {
                    System.Console.WriteLine("Gradient position:");
                    System.Console.WriteLine(item);
                }
                System.Console.WriteLine("\n");
            }


	}
}
