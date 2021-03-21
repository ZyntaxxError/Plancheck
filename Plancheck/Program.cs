////////////////////////////////////////////////////////////////////////////////
//  SRTcheck.cs
//
//  ESAPI v16.1 Script for simple plan parameter checks
//  
////////////////////////////////////////////////////////////////////////////////

using System;
using System.Collections.Generic;
using System.Text.RegularExpressions;
using System.Windows;
using System.Linq;
using VMS.TPS.Common.Model.API;
using VMS.TPS.Common.Model.Types;
using System.IO;
using System.Diagnostics;
using System.CodeDom.Compiler;
using System.Collections;
using System.Text;


/* TODO: 
 * TODO: 
 * TODO Plan sum; foreach plan: check, easier to do in build and wpf...
 * TODO: For dynamic plans; check if verification plan exists only if plan status is planning approved 
 * TODO: z-value for imaging, recommend extended CBCT
 * */

namespace VMS.TPS
{


	public static class Conversion
    {
		public enum Scale
		{
			IEC61217,
			VarianIEC,
			Eclipse,
			Dicom
		}

		public static double ScaleConv(double val, Scale fromScale, Scale toScale, Machine.Axis axis)
        {
            switch (axis)
            {
                case Machine.Axis.Gantry:
                    switch (fromScale)
                    {
                        case Scale.IEC61217:
                            switch (toScale)
                            {
                                case Scale.VarianIEC:
                                    if (val >= 0 && val <= 180 )
                                    {
										val += 180;
                                    }
                                    else
                                    {
										val -= 180;
                                    }
                                    break;
                                case Scale.Eclipse:
                                    break;
                                case Scale.Dicom:
                                    break;
                                default:
                                    break;
                            }
                            break;
                        case Scale.VarianIEC:
                            break;
                        case Scale.Eclipse:
                            break;
                        case Scale.Dicom:
                            break;
                        default:
                            break;
                    }

                    break;
                case Machine.Axis.Coll:
                    break;
                case Machine.Axis.CouchVrt:
                    break;
                case Machine.Axis.CouchLng:
                    break;
                case Machine.Axis.CouchLat:
                    break;
                case Machine.Axis.CouchRtn:
                    break;
                case Machine.Axis.CouchPit:
                    break;
                case Machine.Axis.CouchRol:
                    break;
                case Machine.Axis.y1:
                    break;
                case Machine.Axis.y2:
                    break;
                case Machine.Axis.x1:
                    break;
                case Machine.Axis.x2:
                    break;
                default:
                    break;
            }

            return val;
        }
	}


	// speeds are mean values from trajectory logs, parameters vary somewhat depending on acceleration or deceleration, gravity etc but should be close enough
	// to calculate a fair estimation of beam-on time.
	public static class Machine
	{
		public const int GantryMaxSpeed = 6;  // deg/s
		public const int CollMaxSpeed = 60;  // TODO: TBD
		public const float JawYMaxSpeed = 22.5f;
		public const float JawXMaxSpeed = 22.5f;
		public const int GantryMaxAcc = 16;// deg/s^2 
		public const int CollMaxAcc = 6;  // TODO: TBD
		public const int JawYMaxAcc = 60; // mm/s^2
		public const int JawXMaxAcc = 160; // mm/s^2

		public enum Axis
		{
			Gantry,
			Coll,
			CouchVrt,
			CouchLng,
			CouchLat,
			CouchRtn,
			CouchPit,
			CouchRol,
			Y,
			y1,
			y2,
			X,
			x1,
			x2
		}



		public static readonly Dictionary<string, int> DoseRates = new Dictionary<string, int>
			{
				{ "6X", 600 },
				{ "15X", 600 },
				{ "6XFFF", 1400 },
				{ "10XFFF", 2400 }
			};
	}


		public static class BeamExtraInfo  // TODO: check if this class is better as static methods only, extension class to beam? naming?
	{
		/*
		public BeamExtraInfo(Beam beam)
		{
			Id = beam.Id;
			EBOT = GetEstimatedBeamOnTime(beam);
		}

		public string Id { get; private set; }
		public double EBOT { get; private set; }

		*/

		/// <summary>
		/// Calculates an estimated beam on time from control points in beam, works for dynamic as well as static fields. TODO: EDW unhandled, underestimates 
		/// time for wedges, especially large fields (Y) with small number of MU, need to get the STT-tables and calculate this separately 
		/// </summary>
		/// <param name="beam"></param>
		/// <returns></returns>
		public static double GetEstimatedBeamOnTime(Beam beam)
		{
			double timeOffset = 0.3; // s, added time for startup beam stabilisation, empirical estimation from trajectory logs
			double time = 0;
			int maxDoseRate = beam.DoseRate / 60; // MU/s
			double totalMU = beam.Meterset.Value;
			double deltaMU;

			int nrOfJaws = 4;

			int nrCP = beam.ControlPoints.Count();
			double[] cpTime = new double[nrCP];
			double[] deltaGantry = new double[nrCP];
			double[] gantrySpeed = new double[nrCP];
			double[] doseRate = new double[nrCP];
			double[,] deltaJaw = new double[nrCP, nrOfJaws];
			double[,] jawSpeed = new double[nrCP, nrOfJaws];
			double[] timeJaw = new double[nrOfJaws];

			for (int i = 0; i < nrOfJaws; i++)
			{
				jawSpeed[0, i] = 0;
			}

			double timeGantry;
			gantrySpeed[0] = 0;


			// debug list of control points, calculated time and axis speeds to compare with trajectory logs 
			//TODO: check sign of speeds due to difference in scale (Machine scale vs controlpoints (IEC 61217?)
			StringBuilder controlPointList = new StringBuilder();
			controlPointList.AppendLine("CP\tt\tG_S\tX1_S\tX2_S\tY1_S\tY2_S\tDR");


			// Assumptions: dose rate modulation is instant, MLC speed and acceleration not an issue (might be an issue for banks or for larger MLC leaves, unknown).
			// Time needed from one cp to the next is determined by the parameter that requires the most time, i.e calculate the minimum time needed for each axis

			for (int i = 1; i < beam.ControlPoints.Count(); i++)
			{
				// time needed to deliver the delta MU with specified dose rate
				deltaMU = totalMU * (beam.ControlPoints[i].MetersetWeight - beam.ControlPoints[i - 1].MetersetWeight);
				cpTime[i] = deltaMU / maxDoseRate; // temp assign time for control point

				// minimum time needed for gantry to move between control points at max gantry speed and max acceleration
				deltaGantry[i] = DeltaAngle(beam.ControlPoints[i - 1].GantryAngle, beam.ControlPoints[i].GantryAngle, Conversion.Scale.IEC61217, Machine.Axis.Gantry);
				timeGantry = GetMinTravelTime(deltaGantry[i], gantrySpeed[i - 1], Machine.Axis.Gantry);

				if (timeGantry > cpTime[i])
				{
					cpTime[i] = timeGantry;
				}

				// jaw move time. deltajaw has sign, need to consider sign for speed change
				deltaJaw[i, 0] = beam.ControlPoints[i].JawPositions.X1 - beam.ControlPoints[i - 1].JawPositions.X1;
				deltaJaw[i, 1] = beam.ControlPoints[i].JawPositions.X2 - beam.ControlPoints[i - 1].JawPositions.X2;
				deltaJaw[i, 2] = beam.ControlPoints[i].JawPositions.Y1 - beam.ControlPoints[i - 1].JawPositions.Y1;
				deltaJaw[i, 3] = beam.ControlPoints[i].JawPositions.Y2 - beam.ControlPoints[i - 1].JawPositions.Y2;

				for (int j = 0; j < nrOfJaws; j++)
				{
					if (j < 2)
					{
						timeJaw[j] = GetMinTravelTime(deltaJaw[i, j], jawSpeed[i - 1, j], Machine.Axis.X);     // s
					}
					else
					{
						timeJaw[j] = GetMinTravelTime(deltaJaw[i, j], jawSpeed[i - 1, j], Machine.Axis.Y);     // s
					}
				}

				if (timeJaw.Max() > cpTime[i])
				{
					cpTime[i] = timeJaw.Max();
				}

				// need to recalculate gantry and jaw speed from the resulting control point time to use for next control point
				// TODO: this might be completely off... problably need to study how the control system acts between control points

				gantrySpeed[i] = deltaGantry[i] / cpTime[i];

				for (int j = 0; j < nrOfJaws; j++)
				{
					jawSpeed[i, j] = deltaJaw[i, j] / cpTime[i]; // mm/s
				}

				time += cpTime[i];

				//debug list of control points for comparison with trajectory logs, perhaps get timeJaw?
				controlPointList.AppendLine(i + "\t" + cpTime[i].ToString("0.00") + "\t" + "\tGspeed:" + gantrySpeed[i].ToString("0.0"));
				for (int j = 0; j < nrOfJaws; j++)
				 {
					 controlPointList.Append((jawSpeed[i, j]/10).ToString("0.0") + "\t"); // Convert to cm/s to compare with log files
				}
				doseRate[i] = deltaMU * 60 / cpTime[i];     // calculation of doserate to compare with log files MU/min
				controlPointList.Append("\t" + doseRate[i]);
			}

			string cpInfo = controlPointList.ToString();

			MessageBox.Show(cpInfo);

			return time + timeOffset;
		}


		/// <summary>
		/// Calculates minimum time required to move an axis a distance s. 
		/// </summary>
		/// <param name="s"></param>
		/// <param name="v0"></param>
		/// <param name="axis"></param>
		/// <returns></returns>
		private static double GetMinTravelTime(double s, double v0, Machine.Axis axis)
		{
			double vMax = 1, a = 1, t;
			switch (axis)
			{
				case Machine.Axis.Gantry:
					vMax = Machine.GantryMaxSpeed;
					a = Machine.GantryMaxAcc;
					break;
				case Machine.Axis.Y:
					vMax = Machine.JawYMaxSpeed;
					a = Machine.JawYMaxAcc;
					break;
				case Machine.Axis.X:
					vMax = Machine.JawXMaxSpeed;
					a = Machine.JawXMaxAcc;
					break;
			}

			//The sign of s, vmax and acceleration should always be the same
			int movementDirection = (int)(s / Math.Abs(s));

			vMax *= movementDirection;
			a *= movementDirection;

			// need to solve a second degree equation to get time
			double rot = Math.Sqrt((v0 / a) * (v0 / a) + 2 * s / a);
			t = -v0 / a;
			if (t - rot < 0)
			{
				t += rot;
			}
			else
			{
				t -= rot;
			}

			// if the calculated time is longer than the time it takes to reach vmax, i.e if resulting v > vmax
			if (a * t + v0 > vMax)
			{
				double t1 = (vMax - v0) / a; // time it takes to reach vmax
				double s1 = (a * t1 * t1) / 2 + v0 * t1; // distance traveled when vmax is reached
				t = t1 + (s - s1) / vMax; // add time to travel remaining distance (s-s1) at constant velocity vmax
			}
			return t;
		}


		/// <summary>
		/// Calculates the absolut difference in degrees between two angles considering the scale and axis.
		/// Example: difference between gantry angle 170 and 190 (IEC) should be 140 due to gantry movement restictions. 
		/// TODO: Extended range for gantry angles unhandled, no flag for this in esapi. 
		/// </summary>
		/// <param name="angleStart"></param>
		/// <param name="angleEnd"></param>
		/// <returns></returns>
		private static double DeltaAngle(Double angleStart, Double angleEnd, Conversion.Scale scale, Machine.Axis axis)
		{
			//Convert to Varian machine scale before calculating delta angle 
			angleStart = Conversion.ScaleConv(angleStart, scale, Conversion.Scale.VarianIEC, axis);
			angleEnd = Conversion.ScaleConv(angleEnd, scale, Conversion.Scale.VarianIEC, axis);

			double dAngle = Math.Abs(angleEnd - angleStart);
			
			return dAngle;
		}
	}


	class Script
	{
		public Script()
		{
		}
		// plan categories
		enum PlanCat
		{
			Unknown,
			Electron,
			SBRT,
			CRT,
			IMRT,
			VMAT,
			DynArc,
			TBI,
			TMI,
			NoPlan
		}


		public void Execute(ScriptContext context)
		{
			if ((context.PlanSumsInScope.FirstOrDefault() == null && context.PlanSetup == null) || context.StructureSet == null)
			{
                if (context.StructureSet == null)
                {
					MessageBox.Show("Please select a plan in the active context window.");
				}
                else
                {
					StructureSet ss = context.StructureSet;
					string messageTitle = "Quick check on structure set" + ss.Id;
					string message = CheckStructureSet(ss) +
					CheckCouchStructure(ss);
					MessageBox.Show(message, messageTitle);
				}
			}
			else
			{
                    PlanSetup plan = context.PlanSetup;
                    PlanCat planCat = PlanCat.Unknown;
                    CheckPlanCategory(plan, ref planCat);
                    PlanSum psum = context.PlanSumsInScope.FirstOrDefault();

                    if (psum != null && planCat == PlanCat.TMI)
                    {
                        CheckPlanSum(psum);
                    }
                    else
                    {

                    List<string> relevantDocuments = new List<string>();
					Course course = context.Course;
					string courseIntent = context.Course.Intent;
					string courseId = context.Course.Id;
					string messageTitle = "Quick check on " + courseId + " " + plan.Id;
					string message = string.Empty;

					message +=
						"Assumed plan category: " + planCat + "\n\n" +
					CheckCourseIntent(courseIntent, plan, planCat) + "\n" +
					CheckPlanNamingConvention(plan) +
					CheckClinProt(plan) +
					CheckTargetVolumeID(plan, plan.StructureSet) +
					CheckCtvPtvPlanNumbering(plan, planCat) +
					CheckForBolus(plan, ref relevantDocuments) +
					CheckFieldNamingConvention(course, plan) +
					CheckStructureSet(plan, planCat);


					if (planCat == PlanCat.Electron)
                    {
						message += CheckElectronPlan(plan);
					}
                    else
                    {
						message +=
						CheckSetupField(plan) + "\n" +
						//"Treatment fields \n" +
						//"ID \t Energy \t Tech. \t Drate \t MLC \n" +
						CheckFieldRules(plan, planCat) +
						CheckCouchStructure(plan.StructureSet) +
						CheckStructureSet(plan, planCat);
						//DeltaShiftFromOrigin(plan);
						//SimpleCollisionCheck(plan);
					}


					if (message.Length == 0)
					{
						message = "No comments...";
					}

					MessageBox.Show(message, messageTitle);

				}
			}
		}

        private void CheckPlanSum(PlanSum psum)
        {// Only for TMI...
			string sumPlans = string.Empty;

			//List<PlanSetup> sumPlanTreatOrderHFS = psum.PlanSetups.OrderBy(p => p.TreatmentOrientation).ThenBy(p => p.Beams.Where(b => b.IsSetupField).Count()).ThenBy(p => p.Id).ToList();
			//List<PlanSetup> sumPlanTreatmentOrder = psum.PlanSetups.OrderBy(p => p.Id).ToList();
			List<PlanSetup> sumPlanTreatOrderHFS = psum.PlanSetups.Where(p => p.TreatmentOrientation == PatientOrientation.HeadFirstSupine).OrderBy(p => p.Beams.Where(b => b.IsSetupField).Count()).ThenBy(p => p.Id).ToList();
			List<PlanSetup> sumPlanTreatOrderFFS = psum.PlanSetups.Where(p => p.TreatmentOrientation == PatientOrientation.FeetFirstSupine).OrderBy(p => p.Beams.Where(b => b.IsSetupField).Count()).ThenBy(p => p.Id).ToList();

            if (sumPlanTreatOrderHFS.Count() >= 1)
            {
				sumPlans += "HFS" + "\n\n" + sumPlanTreatOrderHFS[0].Id + "\t\n" + DeltaShiftFromOrigin(sumPlanTreatOrderHFS[0]);

				for (int i = 1; i < sumPlanTreatOrderHFS.Count(); i++)
				{
					sumPlans += sumPlanTreatOrderHFS[i].Id + "\t\n" + DeltaShiftFromPlanToPlan(sumPlanTreatOrderHFS[i - 1], sumPlanTreatOrderHFS[i]);
				}
			}

            if (sumPlanTreatOrderFFS.Count() >= 1)
            {
				sumPlans += "\n\nFFS" + "\n\n" + sumPlanTreatOrderFFS[0].Id + "\t\n" + DeltaShiftFromOrigin(sumPlanTreatOrderFFS[0]);

				for (int i = 1; i < sumPlanTreatOrderFFS.Count(); i++)
				{
					sumPlans += sumPlanTreatOrderFFS[i].Id + "\t\n" + DeltaShiftFromPlanToPlan(sumPlanTreatOrderFFS[i - 1], sumPlanTreatOrderFFS[i]);
				}
			}
			MessageBox.Show(sumPlans);
		}


        #region Plan Category determination


        /// <summary>
        /// Basic check for which plan category the active plan falls into
        /// </summary>
        /// <param name="plan"></param>
        /// <param name="planCat"></param>
        private void CheckPlanCategory(PlanSetup plan, ref PlanCat planCat)
		{
			
			if (IsPlanElectron(plan))
            {
				planCat = PlanCat.Electron;
			}
			else if (IsPlanSRT(plan))
			{
				planCat = PlanCat.SBRT;
			}
			else if (plan.Course.Intent.Contains("TBI"))
            {
				planCat = GetTechniqueTBI(plan);
			}
		}

        private PlanCat GetTechniqueTBI(PlanSetup plan)
        {
			string gd = plan.Beams.Where(b => !b.IsSetupField).First().GantryDirection.ToString();

            if (gd.Contains("None"))
			{
				return PlanCat.TBI;	// Oldschool TBI
            }
            else
            {
				return PlanCat.TMI;
            }
		}

        public bool IsPlanElectron(PlanSetup plan)
		{
			return plan.Beams.Where(b => !b.IsSetupField).Where(b => b.EnergyModeDisplayName.Contains("E")).Any();
		}

		public bool IsPlanSRT(PlanSetup plan)
		{
			bool fractionationSRT = ((plan.NumberOfFractions > 2 && plan.DosePerFraction.Dose >= 8.0 && plan.TotalDose.Dose >= 40.0) || (plan.NumberOfFractions >= 5 && plan.DosePerFraction.Dose >= 6 && plan.TotalDose.Dose >= 35.0));
			bool cIntentSRT = plan.Course.Intent.Contains("SRT");
			return (fractionationSRT || cIntentSRT);
		}


		#endregion Plan Category determination

		#region Electron Plan Checks

		private string CheckElectronPlan(PlanSetup plan)
		{
			//TODO: Check valid SSD
			string cResults = string.Empty;
			foreach (var beam in plan.Beams.Where(b => !b.IsSetupField).OrderBy(b => b.Id)) 
			{
				cResults += beam.Technique.Id + "\n" + beam.Applicator.Id + "\n";
				// Check number of blocks connected to beam
                if (beam.Blocks.Count() == 1)
                {
					Block eBlock = beam.Blocks.First();
					cResults += eBlock.Tray.Id + "\n" + // CustomFFDA, FFDA(A10+)
						eBlock.AddOnMaterial.Id + "\n" + //Elektonblock
						eBlock.IsDiverging + "\n" +
						eBlock.Id + "\n" +
						eBlock.Type + "\n" + //BlockType.APERTURE + "\n";
						CheckEBlockRules(beam) +
						CheckEBlockOutLine(beam, eBlock);
				}
                else if (beam.Blocks.Count() > 1)
                {
					cResults += "WTH!! why would u need more than ONE block?!\n";
				} 
				else
                {
					cResults += "Note: WTF acc. to Centuri, still need a block??????????????No custom block will give the nominal square field size of the applicator with the standard insert + " +
						"Can not automaticly check the Tray Id under Accessories -> Slot 3.\n";
				}
			}
			return cResults;
		}

        private string CheckEBlockRules(Beam beam)
        {
			string cResults = string.Empty;
			Block eBlock = beam.Blocks.First();
			if (eBlock.Type != BlockType.APERTURE)
			{
				cResults += "* Electron block type should always be Aperture, not Shielding!\n";
			}
			if (eBlock.IsDiverging)
			{
				cResults += "* Electron block type should generally not be diverging.\n";
			}
			if (beam.Applicator.Id.Contains("A06"))
			{
				cResults += "* Please choose appl. A10 instead if custom block is needed, or remove the block if nominal field size is desired.\n";
			}
			else if (!eBlock.Tray.Id.Contains("CustomFFDA"))
			{
				// THIS IS APPARANTLY WRONG!?
				cResults += "** Wrong Block Tray selected! See guidelines in Centuri.\n";
			}
			return cResults;
		}

        private string CheckEBlockOutLine(Beam beam, Block eBlock)
        {
			string cResults = string.Empty;
			var blockPoints = eBlock.Outline;
			int[] standardEBlockDiameters = { 4, 5, 6, 7, 8, 9, 10 };
			double radii = 0;
            double minRadi = 1000;
            double maxRadi = 0;
			// TODO: I assume this whole thing can be made into a one liner with linq...
            foreach (var p in blockPoints)
            {
                for (int i = 0; i < p.Count(); i++)
                {
					Vector pointVector = (Vector)p[i];
					radii = pointVector.Length;
                    if (radii > maxRadi)
                    {
						maxRadi = radii;
                    }
                    else if (radii < minRadi)
                    {
						minRadi = radii;
                    }
                }
			}
            if (Math.Round(maxRadi, 0) == Math.Round(minRadi, 0) && standardEBlockDiameters.Contains((int)(2 * Math.Round(maxRadi, 0)/10)))
			{
				int stdBlockDiam = (int)(2 * Math.Round(maxRadi, 0) / 10);

				if (beam.Applicator.Id.Equals("A10"))
                {
					// TODO: check the naming of the block ID, which should reflect the diameter of the block
					cResults += "Standard block with diameter " + stdBlockDiam.ToString("0") + " cm\n";
					cResults += CheckStandardEBlockNamingConvention(eBlock.Id, stdBlockDiam);
				}
                else
                {
					cResults += "* The block aperture diameter (" + ((int)(2 * Math.Round(maxRadi, 0)) / 10).ToString("0") + " cm)" +
						" equals a standard block, please change the applicator to A10 if you want to use a standard block!\n";
				}
			}
			return cResults;
		}


		private string CheckStandardEBlockNamingConvention(string blockID, int apertureDiameter)
		{
			string cResults = string.Empty;
			int blockIdDiam = 0;
			string wrongNaming = "* Block ID should include the aperture diameter for standard blocks \n( i.e. \"Diameter " + apertureDiameter.ToString("0") + " cm\" ).";

			var blockApertureDiamCm = new Regex(@"(\d{1,2})\s?(cm)", RegexOptions.IgnoreCase);
			var blockApertureDiamMm = new Regex(@"(\d{1,3})\s?(mm)", RegexOptions.IgnoreCase);

			if (blockApertureDiamCm.IsMatch(blockID))
			{
				Group g = blockApertureDiamCm.Match(blockID).Groups[1];
				if (int.TryParse(g.Value, out blockIdDiam))
				{
					cResults += "";
				}
			}
			else if (blockApertureDiamMm.IsMatch(blockID))
			{
				Group g = blockApertureDiamMm.Match(blockID).Groups[1];
				if (int.TryParse(g.Value, out blockIdDiam))
				{
					blockIdDiam = (int)(blockIdDiam / 10); 
				}
            }
            else
            {
				cResults += wrongNaming;
			}
            if (apertureDiameter != blockIdDiam)
            {
				cResults += wrongNaming;
			}

			return cResults;
		}





		#endregion Electron Plan Checks



		private string CheckCtvPtvPlanNumbering(PlanSetup plan, PlanCat planCat)
		{
			string cResults = string.Empty;
			StructureSet ss = plan.StructureSet;
			int number;

			// check that CTV, PTV and Plan numbering is consistent, only nesessary if more than one PTV
			// PTV from plans target volume, ctv from PTV number and check that mass center is within respective PTV boundaries
			// Plan numbering not necessarily consistent with PTV if multiple structure sets used in same course 
			if (!string.IsNullOrEmpty(plan.TargetVolumeID) && planCat == PlanCat.SBRT)
			{
				// Search for structure in structure set with same id as target volume and checks if type is PTV, defaults to null if criteria not met
				Structure ptv = ss.Structures.Where(s => s.Id == plan.TargetVolumeID).Where(s => s.DicomType == "PTV").SingleOrDefault();
				// Check if more than one PTV exists
				int nrOfPTV = ss.Structures.Where(s => s.Id.Length > 4).Where(s => s.Id.Substring(0, 3) == "PTV").Where(s => s.DicomType == "PTV").Count();
				if (ptv != null && Int32.TryParse(ptv.Id.Substring(3, 1), out number) && nrOfPTV > 1)
				{
					//Structure ctv = ss.Structures.Where(s => s.Id.Substring(0, 4) == ("CTV" + ptv.Id.Substring(3, 1))).Where(s => s.DicomType == "CTV").SingleOrDefault();
					if (plan.Id.Substring(1, 1) != ptv.Id.Substring(3, 1))
					{
						cResults += "* Plan ID number not consistent with PTV ID number, check if this is ok";
					}
					//cResults = "Plan,PTV and CTV \t" + plan.Id + "\t" + ptv.Id + "\t" + ctv.Id + "\t" + nrOfPTV + "\n";
				}
			}
			return cResults + "\n";
		}


		#region delta shift

		private string DeltaShiftFromPlanToPlan(PlanSetup fromPlan, PlanSetup toPlan)
		{
			VVector deltaShift = DeltaShiftIncm(fromPlan, fromPlan.Beams.First().IsocenterPosition, toPlan.Beams.First().IsocenterPosition);


			string delta = " Delta(Vrt,Lng,Lat)[cm]: \t" + deltaShift.y.ToString("0.00");
			delta += "\t" + deltaShift.z.ToString("0.00");
			delta += "\t" + deltaShift.x.ToString("0.00")+ "\n";
			return delta;
		}





		private string DeltaShiftFromOrigin(PlanSetup plan)
        {
			VVector deltaShift = DeltaShiftIncm(plan, plan.StructureSet.Image.UserOrigin, plan.Beams.First().IsocenterPosition);


			string delta = "\niso-Vrt: \t" + deltaShift.y.ToString("0.00") + " cm\n";
			delta += "iso-Lng: \t" + deltaShift.z.ToString("0.00") + " cm\n";
			delta += "iso-Lat: \t" + deltaShift.x.ToString("0.00") + " cm\n\n";
			return delta;
		}

        private VVector DeltaShiftIncm(PlanSetup plan, VVector dicomOriginalPosition, VVector dicomFinalPosition)
        {
			Image image = plan.StructureSet.Image;
			VVector eclipseOriginalPosition = image.DicomToUser(dicomOriginalPosition, plan);
			VVector eclipseFinalPosition = image.DicomToUser(dicomFinalPosition, plan);
			VVector deltaShift = (eclipseOriginalPosition - eclipseFinalPosition) / 10;

			// Have to change sign of Vrt before returning due to difference in coordinate system in Eclipse and the machine
			deltaShift.y *= -1;

			return deltaShift;
		}

		#endregion

		#region Bolus // TODO: thickness in ID probably deprecated in 16.x...


		// ********* Kontroll av bolus, om kopplat till alla fält, förväntat HU-värde och namngivning  *********
		// kollar enbart bolus kopplat till något fält
		private string CheckForBolus(PlanSetup plan, ref List<string> relevantDocuments)
        {
			string cResults = string.Empty;
			List<string> bolusList = new List<string>();
			List<int> bolusListnr = new List<int>();
			List<double> bolusListHU = new List<double>();
			int noOfBoluses = 0;
			int nrOfBeams = plan.Beams.Count();
			// iterate through all beams in plan, and then all boluses in beams to get all boluses used
            foreach (var beam in plan.Beams)
            {
                foreach (var bolus in beam.Boluses)
                {
					noOfBoluses += beam.Boluses.Count();
                    if (bolusList.IndexOf(bolus.Id) < 0)
                    {
						bolusList.Add(bolus.Id);
						bolusListnr.Add(1);
						bolusListHU.Add(bolus.MaterialCTValue);
					}
                    else
                    {
						bolusListnr[bolusList.IndexOf(bolus.Id)] += 1;
					}
                }
            }
            if (noOfBoluses >= 1)
            {
				relevantDocuments.Add("https://www.red-gate.com/si");
			}
            for (int i = 0; i < bolusList.Count(); i++)
            {
                if (bolusListnr[i] != nrOfBeams)
                {
					cResults += "* Bolus \"" + bolusList[i] + "\" not linked to all fields, check if this is intended. \n\n";
				}
                if (bolusListHU[i] != 0)
                {
					cResults += "* Bolus \"" + bolusList[i] + "\" not assigned default HU (assigned " + bolusListHU[i].ToString("0") + " HU instead of expected/default 0 HU), check if this is intended. \n\n";
				}
				cResults += CheckBolusNamingConvention(bolusList[i]);
            }
			return cResults;
        }

        private string CheckBolusNamingConvention(string bolusID)
		{ 
			string cResults = string.Empty;
			double bolusThick = 0;

			var bolusThicknessCm = new Regex(@"(\d{1}.?\d?)\s?(cm)", RegexOptions.IgnoreCase);
			var bolusThicknessMm = new Regex(@"(\d{1,2})\s?(mm)", RegexOptions.IgnoreCase);
			var bolusCTpattern = new Regex(@"ct", RegexOptions.IgnoreCase);
            if (bolusThicknessCm.IsMatch(bolusID))
            {
				Group g = bolusThicknessCm.Match(bolusID).Groups[1];
				if (Double.TryParse(g.Value.Replace(",", ".") , out bolusThick))
                {
					cResults += CheckBolusThickness((int)(bolusThick * 10));
                }
            }
            else if (bolusThicknessMm.IsMatch(bolusID))
            {
				Group g = bolusThicknessMm.Match(bolusID).Groups[1];
				if (Double.TryParse(g.Value, out bolusThick))
				{
					cResults += CheckBolusThickness((int)bolusThick);
				}
            }
            else if (!bolusCTpattern.IsMatch(bolusID))
            {
				cResults += "* Check naming convention of bolus ID: \"" + bolusID + "\"\n\n";
			}

            return cResults;
        }

        private string CheckBolusThickness(int mmBolus)
        {
			string cResults = string.Empty;
			int[] availableThicknesses = { 3, 5, 8, 10, 13, 15, 18, 20 };
            if (!availableThicknesses.Contains(mmBolus))
            {
				cResults = "* Physical thickness of available boluses are 3, 5 and 10 mm and a combination thereof. \n\n";
			}
			return cResults;
        }

        #endregion Bolus

        #region Structure set

        // TMI: order plans from head to toe (by isopos in dicom, reverse), check ID numbering
        // search for mask base plate and half spheres in image


        private string CheckStructureSet(PlanSetup plan, PlanCat planCategory)
		{
			StructureSet ss = plan.StructureSet;
			string cResults = string.Empty;
			string cAssignedHU = string.Empty;
			try
            {
				foreach (var structure in ss.Structures.Where(s => !s.Id.Contains("Couch")))
				{
					cAssignedHU += CheckForAssignedHU(structure);
				}
				if (planCategory == PlanCat.SBRT)
				{
					cResults += CheckStructureSetSBRT(plan);
				}
			}
            catch (Exception e)
            {
				cResults += e;
			}
			if (!string.IsNullOrEmpty(cAssignedHU))
			{
				cResults += "\nStructures with assigned HU: \n" + cAssignedHU + "\n\n";
			}

			return cResults;
		}

		private string CheckStructureSet(StructureSet ss)
		{
			string cResults = string.Empty;
			string cAssignedHU = string.Empty;
			string cSeparateParts = string.Empty;
			int numberOfSeparateParts;
			try
			{
				foreach (var structure in ss.Structures.Where(s => !s.Id.Contains("Couch")).Where(s => !s.Id.Contains("Bones")))
				{
					cAssignedHU += CheckForAssignedHU(structure);
                    if (!structure.IsEmpty && structure.HasSegment)
                    {
						numberOfSeparateParts = structure.GetNumberOfSeparateParts();
						if (numberOfSeparateParts > 1)
						{
							cSeparateParts += string.Format("Structure \t{0} has \t{1} separate parts.\n", structure.Id, numberOfSeparateParts);
						}
					}
				}
			}
			catch (Exception e)
			{
				cResults += e;
			}
            if (!string.IsNullOrEmpty(cAssignedHU))
            {
				cResults += "Structures with assigned HU: \n\n" + cAssignedHU + "\n\n";
			}
			if (!string.IsNullOrEmpty(cSeparateParts))
			{
				cResults += "Structures with separate parts: \n\n" + cSeparateParts + "\n\n";
			}


			// attemt to roughly estimate CTV to PTV margin
			// Search for structure in structure set with same id as target volume and checks if type is PTV, defaults to null if criteria not met
			List<Structure> ptvList = ss.Structures.Where(s => s.Id.Substring(0,3) == "PTV").Where(s => s.DicomType == "PTV").ToList();
			// Check if more than one PTV exists
			if (ptvList != null)
			{
				cResults += "Estimation of CTV to PTV margin by comparing outer bounds of respective structure:\n";
				foreach (var ptv in ptvList)
				{
                    try
                    {
						Structure ctv = ss.Structures.Where(s => ptv.Id.Contains(s.Id.Substring(1, s.Id.Length - 1))).Where(s => s.DicomType == "CTV").SingleOrDefault();
						//Structure ctv = ss.Structures.Where(s => ptv.Id.Contains(s.Id.Substring(1, s.Id.Length - 1))).SingleOrDefault();
						bool correctCTV = true;
						if (ctv == null)
						{
							ctv = ss.Structures.Where(s => s.Id.Substring(0, 3) == "CTV").Where(s => (s.CenterPoint - ptv.CenterPoint).Length < 5).Where(s => ptv.MeshGeometry.Bounds.Contains(s.MeshGeometry.Bounds)).SingleOrDefault();
						}
						else
						{
							// if identical ID, check that PTV actually contains CTV 
							if ((ctv.CenterPoint - ptv.CenterPoint).Length > 5 || !ptv.MeshGeometry.Bounds.Contains(ctv.MeshGeometry.Bounds))
							{
								//ctv = null; Won't allow this
								correctCTV = false;
							}
						}
						//Structure ctv = ss.Structures.Where(s => s.Id.Substring(1, s.Id.Length - 1) == ptv.Id.Substring(1, ptv.Id.Length - 1)).SingleOrDefault();
						// unfortunately are some ctv in clinical protocol dicomType avoid...
						// strongle depends of exactly identical naming, can perhaps compare mass centrum or bounds contain?
						if (ctv != null && correctCTV)
						{
							double deltaX = Math.Round(ptv.MeshGeometry.Bounds.SizeX + 0.5 - ctv.MeshGeometry.Bounds.SizeX) / 2;
							double deltaY = Math.Round(ptv.MeshGeometry.Bounds.SizeY + 0.5 - ctv.MeshGeometry.Bounds.SizeY) / 2;
							double deltaZ = Math.Round(ptv.MeshGeometry.Bounds.SizeZ - ctv.MeshGeometry.Bounds.SizeZ) / 2;
							cResults += ctv.Id + " -> " + ptv.Id + ":\t\tx: " + deltaX.ToString("0") + "\ty: " + deltaY.ToString("0") + "\tz: " + deltaZ.ToString("0") + "\n";
						}
						else
						{
							cResults += ctv.Id + " -> " + ptv.Id + ":\t can not identify corresponding CTV.\n";
						}
					}
                    catch (Exception e)
                    {

						cResults += ptv.Id + " gave following error: \t" + e;

					}

				}
			}
			return cResults;
		}


		private string CheckForAssignedHU(Structure structure)
        {
			string cResults = string.Empty;
			double huValue;
			if (structure.GetAssignedHU(out huValue))
			{
				cResults = string.Format("{0} has an assigned CT value of \t{1} HU\n", structure.Id, huValue);
			}
			return cResults;
		}

        private string CheckStructureSetSBRT(PlanSetup plan)
        {
			StructureSet ss = plan.StructureSet;
			string cResults = string.Empty;

			Image image = ss.Image;
            if (image.ZRes >= 2.5)
            {
				cResults += "* Image slice thickness is " + image.ZRes.ToString("0.0") + " mm. Check if this is ok for this SBRT case.\n";

			}
			Structure skin = ss.Structures.Where(s => s.Id == "Skin").SingleOrDefault();
			if (skin == null || skin.IsEmpty)
			{
				cResults = "* missing structure Skin \n";
            }
            else
            {


				// check that body is larger than Skin, at least in iso pos (z)
				// Could use bounds but that assumes that skin is contured for whole body 
				// Instead use profile
				// startpoint; isocenter long and lat, vrt bottom of bounding box for total body

				// If vacum bag, could be thin in bottom, better get a profile also in lateral direction

				Structure body = ss.Structures.Where(s => s.Id == "BODY").SingleOrDefault();  // TODO: better to take dicom type
				double bodyVrtMin = body.MeshGeometry.Bounds.Y + body.MeshGeometry.Bounds.SizeY;

				VVector endPoint = plan.Beams.FirstOrDefault().IsocenterPosition;
				VVector startPoint = endPoint;
				startPoint.y = bodyVrtMin;
				

				var bodyEntrance = body.GetSegmentProfile(startPoint, endPoint, new BitArray(100)).Where(x => x.Value == true).First();
				var skinEntrance = skin.GetSegmentProfile(startPoint, endPoint, new BitArray(100)).Where(x => x.Value == true).First();

                if (bodyEntrance.Position.y <= skinEntrance.Position.y + 3)
                {
					cResults += "\n* Check if there are any fixation devices that should be included in the BODY structure.\n";
                }
			}
			return cResults;
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



        #endregion



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

		private string CheckCourseIntent(string cIntent, PlanSetup plan, PlanCat planCat)
		{
			string cResults = "";
			if (cIntent.Trim().Length == 0)
			{
				cResults = "** Course intent is empty!";
			}
			else if (planCat == PlanCat.SBRT && !plan.Course.Intent.Contains("SRT"))
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
				Structure target = sSet.Structures.Where(s => s.Id == plan.TargetVolumeID).Where(s => s.DicomType == "PTV" || s.Id.Contains("(fg)")).SingleOrDefault();
				if (target == null)
				{
					cResults = "* Plan target volume should be of type PTV. \n";
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



        #region Set-up field checks




        // ********* 	Kontroll av Setup-fält; namngivning och ej bordsvridning ********* 

        private string CheckSetupField(PlanSetup plan)
		{
			string cResults = "";
			int countSetupfields = 0;
			foreach (var beam in plan.Beams.Where(b => b.IsSetupField))
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
						cResults = cResults + CheckPlanarSetupField(beam, plan) + "\n";    // extra checks for planar setup fields
					}
				}
				else
				{
					cResults = cResults + "** ID should start with Plan-ID" + "\n";
				}
			}
			if (countSetupfields == 0)
			{
				cResults = cResults + "missing!" + "\t \t" + "** Insert Setup field" + "\n";
			}
			// If case there is only one planar setup field: If treatment angle != 0 or 180, error, else reminder to either use catalyst or add a second setup field
			if (plan.Beams.Where(b => b.IsSetupField).Where(b => !b.Id.ToUpper().Contains("CBCT")).Where(b => !b.Id.ToUpper().Contains("BOLUS")).Count() == 1) // CBCT and bolus helper field is ok
			{
                if (plan.Beams.Where(b => b.IsSetupField == false).Where(b => !b.ControlPoints.FirstOrDefault().GantryAngle.Equals(0.0) && !b.ControlPoints.FirstOrDefault().GantryAngle.Equals(180.0)).Count() > 0)
                {
					cResults += "* Only one planar setup field found. Add an orthogonal setup field for kV/kV matching.\n";
                }
                else if(plan.ApprovalStatus == PlanSetupApprovalStatus.PlanningApproved || plan.ApprovalStatus == PlanSetupApprovalStatus.TreatmentApproved)
                {
					cResults += " Only one planar setup field found, remember to add Catalyst for setup (or add another setup field).\n";
				}
			}
			return cResults;
		}

		// ********* 	Kontroll av planara Setup-fält (ej namngivna cbct), namngivning efter gantryvinkel	********* 

		public string CheckPlanarSetupField(Beam beam, PlanSetup plan)
		{
			string cResults = "";
			string trimmedID = beam.Id.Substring(2).Trim();         // start iteration at index 2 (PX index 0 and 1, already checked)
			int gantryAngleInBeamID = 1000;                         // dummy number
			int test = 1000;
			// check for naming according to gantry angle
			for (int i = 1; i < trimmedID.Length + 1; i++)
			{
				if (Int32.TryParse(trimmedID.Substring(0, i), out test))                    //  step up one char at a time and try parsing it to int. If successful assign it to gantryAngleInBeamID
				{
					gantryAngleInBeamID = test;
				}
			}

			// simple collision check only for the major axes
			if (gantryAngleInBeamID != 1000)
			{
				if (gantryAngleInBeamID == Math.Round(beam.ControlPoints.First().GantryAngle))
				{
					cResults += CheckSetupFieldForCollision(beam, plan);
				}
				else
				{
					cResults += "* Check name (gantry angle)!";
				}
			}
			else
			{
                if (beam.Id.ToUpper().Contains("BOLUS"))
                {
					cResults += "OK";
                }
                else
                {
					cResults += "* Id should include gantry angle!";
				}
			}
			return cResults;
		}

        private string CheckSetupFieldForCollision(Beam beam, PlanSetup plan)
        {
			string cResults = string.Empty;
			int latLimitForCollisionCheck = 40;
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
			return cResults;
		}


		#endregion

		// ********* 	Kontroll av diverse fältregler och "best practices"	********* 
		// TODO: better sorting, maybe one general and then divided in categories (enum). Missing case for Static-I and MLC doseDynamic (IMRT)

		private string CheckFieldRules(PlanSetup plan, PlanCat planCat )
		{
			string cResults = string.Empty;
			string remarks = string.Empty;
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
			double beamOnTimeInSec = 0;




			if (planCat == PlanCat.TMI)
			{
				remarks += CheckForCouchValuesTMI(plan);
			}





			foreach (var beam in plan.Beams.Where(b => !b.IsSetupField).OrderBy(b => b.Id))
			{
				//BeamExtraInfo beamExtras = new BeamExtraInfo(beam);

				cResults = cResults + beam.Id + "\t" + beam.EnergyModeDisplayName + "\t" + beam.Technique.Id + "\t" + beam.DoseRate + "\t" + beam.MLCPlanType + "\t";
				cResults += Math.Round(BeamExtraInfo.GetEstimatedBeamOnTime(beam)).ToString("0") + " s\t" + (beam.ControlPoints.Count()) +"\n";
				beamOnTimeInSec += GetEstimatedBeamOnTime(beam);
				


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
			}
			if (Math.Abs(countArcCCW - countArcCW) > 1)
			{
				remarks += "** Check the arc directions! \t";
			}
			return cResults + "\n" + "Estimated total beam-on-time: " + (beamOnTimeInSec/60).ToString("0.0") + " min\n\n" + "TreatTime: > " + Math.Round(GetEstimatedTreatmentTime(plan)/60, 1).ToString("0.0") + remarks;
		}

        private string CheckForCouchValuesTMI(PlanSetup plan)
        {
			string cResults = string.Empty;
			string beamIdWithCouchValues = string.Empty;
            foreach (var beam in plan.Beams)
            {
              if (!beam.ControlPoints[0].TableTopLateralPosition.ToString().Contains("NaN") || !beam.ControlPoints[0].TableTopLongitudinalPosition.ToString().Contains("NaN") || !beam.ControlPoints[0].TableTopVerticalPosition.ToString().Contains("NaN"))
                {
					beamIdWithCouchValues += beam.Id + ", " ;
				}
			}

            if (beamIdWithCouchValues.Length > 0)
            {
				cResults = "* Please remove couch coordinates from the following fields: " + beamIdWithCouchValues.Substring(0, beamIdWithCouchValues.Length - 2) + ".\n\n";
            }
			return cResults;
        }


		/// <summary>
		/// Calculates an estimated beam on time from control points in beam, works for dynamic as well as static fields. TODO: EDW unhandled, underestimates 
		/// time for wedges, espesially large fields (Y) with small number of MU, need to get the STT-tables and calculate this separately 
		/// </summary>
		/// <param name="beam"></param>
		/// <returns></returns>
		private double GetEstimatedBeamOnTime(Beam beam)
		{
			double timeOffset = 0.3; // s, added time for startup beam stabilisation, empirical estimation from trajectory logs
			double time = 0;
			int maxDoseRate = beam.DoseRate / 60; // MU/s
			double totalMU = beam.Meterset.Value;
			double deltaMU;

			

			int nrOfJaws = 4;

			int nrCP = beam.ControlPoints.Count();
			double[] cpTime = new double[nrCP];
			double[] deltaGantry = new double[nrCP];
			double[] gantrySpeed = new double[nrCP];
			double[] jawSpeedX1 = new double[nrCP];
			double[] jawSpeedX2 = new double[nrCP];
			double[] jawSpeedY1 = new double[nrCP];
			double[] jawSpeedY2 = new double[nrCP];
			double[] doseRate = new double[nrCP];

			double[,] deltaJaw = new double[nrCP, nrOfJaws];
			double[,] jawSpeed = new double[nrCP, nrOfJaws];
			double[] timeJaw = new double[nrOfJaws];
			//double[] deltaJawSpeed = new double[nrOfJaws];
            for (int i = 0; i < nrOfJaws; i++)
            {
				jawSpeed[0, i] = 0;
            }
			
			double timeGantry;
			gantrySpeed[0] = 0;
			jawSpeedX1[0] = 0;
			jawSpeedX2[0] = 0;
			jawSpeedY1[0] = 0;
			jawSpeedY2[0] = 0;


			StringBuilder controlPointList = new StringBuilder();

			controlPointList.AppendLine("CP\tt\tdG\tGspeed\tdX1\tdX2\tdY1\tdY2");


			// Assumptions: dose rate modulation is instant. MLC-speed and acceleration not an issue (could be an issue for banks, also for larger MLC leaves, unknown).
			// time needed from one cp to the next is determined by the parameter that requires the most time, i.e calculate the minimum time needed for each axis

			for (int i = 1; i < beam.ControlPoints.Count(); i++)
            {

				// time needed to deliver the delta MU with specified dose rate
				deltaMU = totalMU * (beam.ControlPoints[i].MetersetWeight - beam.ControlPoints[i - 1].MetersetWeight);
				cpTime[i] = deltaMU / maxDoseRate; // temp assign time for control point

				// minimum time needed for gantry to move between control points at max gantry speed including term for acceleration
				deltaGantry[i] = DeltaAngle(beam.ControlPoints[i].GantryAngle,  beam.ControlPoints[i - 1].GantryAngle);
				timeGantry = GetMinTravelTime(deltaGantry[i], gantrySpeed[i - 1], Machine.Axis.Gantry);

				if (timeGantry > cpTime[i]) 
                {
					cpTime[i] = timeGantry;
                }

				// jaw move time. deltajaw has sign, need to consider sign for speed change
				deltaJaw[i, 0] = beam.ControlPoints[i].JawPositions.X1 - beam.ControlPoints[i - 1].JawPositions.X1;
				deltaJaw[i, 1] = beam.ControlPoints[i].JawPositions.X2 - beam.ControlPoints[i - 1].JawPositions.X2;
				deltaJaw[i, 2] = beam.ControlPoints[i].JawPositions.Y1 - beam.ControlPoints[i - 1].JawPositions.Y1;
				deltaJaw[i, 3] = beam.ControlPoints[i].JawPositions.Y2 - beam.ControlPoints[i - 1].JawPositions.Y2;

				for (int j = 0; j < nrOfJaws; j++)
				{
                    if (j<2)
                    {
						timeJaw[j] = GetMinTravelTime(deltaJaw[i, j], jawSpeed[i - 1, j], Machine.Axis.X);     // s
					}
                    else
                    {
						timeJaw[j] = GetMinTravelTime(deltaJaw[i, j], jawSpeed[i - 1, j], Machine.Axis.Y);     // s
					}
				}

                if (timeJaw.Max() > cpTime[i])
                {
					cpTime[i] = timeJaw.Max();
                }


				// need to recalculate gantry and jaw speed from the resulting control point time to use for next control point
				// this might be completely off... problably need to study how the control system acts between control points

				gantrySpeed[i] = deltaGantry[i] / cpTime[i];

                for (int j = 0; j < nrOfJaws; j++)
                {
					jawSpeed[i, j] = deltaJaw[i, j] / cpTime[i];
                }
				doseRate[i] = deltaMU * 60 / cpTime[i];		// calculation of doserate to compare with log files MU/min

				time += cpTime[i];


				controlPointList.AppendLine(i + "\t" + cpTime[i].ToString("0.00") + "\t" + deltaGantry[i].ToString("0.0") + "\tGspeed:" + gantrySpeed[i].ToString("0.0") + "\t" + doseRate[i]);
               /* for (int j = 0; j < nrOfJaws; j++)
                {
					controlPointList.Append(deltaJaw[i, j].ToString("0.0") + "\t");
                }*/
			}

			string cpInfo = controlPointList.ToString();

			MessageBox.Show(cpInfo);

			return time + timeOffset;
		}


		/// <summary>
		/// Calculates minimum time required to move an axis a distance s. 
		/// </summary>
		/// <param name="s"></param>
		/// <param name="v0"></param>
		/// <param name="axis"></param>
		/// <returns></returns>
        private double GetMinTravelTime(double s, double v0, Machine.Axis axis)
        {
			double vMax = 1, a = 1, t;
            switch (axis)
            {
                case Machine.Axis.Gantry:
					vMax = Machine.GantryMaxSpeed;
					a = Machine.GantryMaxAcc;
                    break;
                case Machine.Axis.Y:
					vMax = Machine.JawYMaxSpeed;
					a = Machine.JawYMaxAcc;
					break;
                case Machine.Axis.X:
					vMax = Machine.JawXMaxSpeed;
					a = Machine.JawXMaxAcc;
					break;
            }

			//The sign of s, vmax and acceleration should always be the same

			int movementDirection = (int)(s / Math.Abs(s));

			vMax *= movementDirection;
			a *= movementDirection;

			double rot = Math.Sqrt((v0 / a) * (v0 / a) + 2 * s / a);
			t = -v0 / a;
			if (t-rot < 0)
            {
				t += rot;
            }
            else
            {
				t -= rot;
            }

			// if the calculated time is longer than the time it takes to reach vmax, i.e if resulting v > vmax

            if (a * t + v0 > vMax)
            {
				double t1 = (vMax - v0) / a; // time it takes to reach vmax
				double s1 = (a * t1 * t1) / 2 + v0 * t1; // distance traveled when vmax is reached

				t = t1 + (s-s1) / vMax; // add time to travel remaining distance (s-s1) at constant velocity vmax
            }
			return t;
		}




        /// <summary>
        /// Estimates treatment time (excl. imaging) when using field order according to Beam ID
        /// </summary>
        /// <param name="plan"></param>
        /// <returns></returns>
        private double GetEstimatedTreatmentTime(PlanSetup plan)
        {
			double controlSeqTime;
			double treatTime = 0;
			double maxGantrySpeed = 6;      // deg/s
			double maxJawSpeed = 12.5;      // mm/s
			double manualBeamOnTime = 6;    // s , estimation of time needed to manually push beam on, as opposed to automation

			List<double> mechMovementTime = new List<double>();
			List<ControlPoint> cpFirst = new List<ControlPoint>();
			List<ControlPoint> cpLast = new List<ControlPoint>();
			List<Beam> beamsInOrder = new List<Beam>();
			foreach (var beam in plan.Beams.Where(b => !b.IsSetupField).OrderBy(b => b.Id))
            {
				cpFirst.Add(beam.ControlPoints.First());
				cpLast.Add(beam.ControlPoints.Last());
				beamsInOrder.Add(beam);
				treatTime += GetEstimatedBeamOnTime(beam);
				
            }

			// control sequens and energy change happens parallel with, and independent to, mechanical movements
            for (int i = 1; i < cpFirst.Count(); i++)
            {
				mechMovementTime.Add(DeltaAngle(cpLast[i - 1].GantryAngle, cpFirst[i].GantryAngle)/maxGantrySpeed);
				mechMovementTime.Add(MaxDeltaJaw(cpLast[i - 1], cpFirst[i]) / maxJawSpeed);
                if (beamsInOrder[i-1].EnergyModeDisplayName.Equals(beamsInOrder[i]))
                {
					controlSeqTime = 10;
                }
                else
                {
					controlSeqTime = 20;
                }

                if (controlSeqTime > mechMovementTime.Max())
                {
					treatTime += controlSeqTime;
                }
                else
                {
					treatTime += mechMovementTime.Max();
				}

				mechMovementTime.Clear();

				if (AutomationPrerequisites(plan) == false)
				{
					treatTime += manualBeamOnTime;
				}

			}
				return treatTime;

        }


		private double DeltaAngle(Double angle1, Double angle2)
        {
			double dAngle = Math.Abs(angle2 - angle1);

			if (dAngle > 180)
            {
				return 360 - dAngle;
			}
            else
            {
				return dAngle;
			}
        }



		/// <summary>
		/// gets the maximum distance any jaw travels between two control points in mm.
		/// </summary>
		/// <param name="cp2"></param>
		/// <param name="cp1"></param>
		/// <returns></returns>
		private double MaxDeltaJaw(ControlPoint cp2, ControlPoint cp1)
		{
			double[] deltaJaw = new double[4];

				deltaJaw[0] = Math.Abs(cp2.JawPositions.X1 - cp1.JawPositions.X1);
				deltaJaw[1] = Math.Abs(cp2.JawPositions.X2 - cp1.JawPositions.X2);
				deltaJaw[2] = Math.Abs(cp2.JawPositions.Y1 - cp1.JawPositions.Y1);
				deltaJaw[3] = Math.Abs(cp2.JawPositions.Y2 - cp1.JawPositions.Y2);

			return deltaJaw.Max();
		}




		private bool AutomationPrerequisites(PlanSetup plan) 
		{
			//TODO: Missing check for same accessories for all beams
            if (SingleIsoInPlan(plan) && SingleCouchPosInPlan(plan) && SingleEnergyInPlan(plan))
            {
				return true;
            }
            else
            {
				return false;
            }
		}


		private bool SingleIsoInPlan(PlanSetup plan)
        {
			bool sIso = true;
			Beam firstBeam = plan.Beams.First();
			foreach (var beam in plan.Beams)
			{
				if (beam.IsocenterPosition.x != firstBeam.IsocenterPosition.x || beam.IsocenterPosition.y != firstBeam.IsocenterPosition.y || beam.IsocenterPosition.z != firstBeam.IsocenterPosition.z)
				{
					sIso = false;
					break;
				}
			}
			return sIso;
		}

		private bool SingleEnergyInPlan(PlanSetup plan)
		{
			bool sEnergy = true;
			Beam firstBeam = plan.Beams.Where(b => !b.IsSetupField).OrderBy(b => b.Id).First();
			foreach (var beam in plan.Beams.Where(b => !b.IsSetupField).OrderBy(b => b.Id))
			{
				if (!beam.EnergyModeDisplayName.Equals(firstBeam.EnergyModeDisplayName))
				{
					sEnergy = false;
					break;
				}
			}
			return sEnergy;
		}

		private bool SingleCouchPosInPlan(PlanSetup plan)
		{
			bool sCouch = true;
			Beam firstBeam = plan.Beams.Where(b => !b.IsSetupField).OrderBy(b => b.Id).First();
            if (firstBeam.ControlPoints.First().PatientSupportAngle >= 345 || firstBeam.ControlPoints.First().PatientSupportAngle <= 15)
            {
				foreach (var beam in plan.Beams.Where(b => !b.IsSetupField).OrderBy(b => b.Id))
				{
					if (beam.ControlPoints.First().PatientSupportAngle != firstBeam.ControlPoints.First().PatientSupportAngle)
					{
						sCouch = false;
						break;
					}

				}
            }
            else
            {
				sCouch = false;
			}
			return sCouch;
		}


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

		// ********* 	Kontroll av energi vid Statisk MLC och SBRT	*********

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
				cResults = "* Collimator angle for DynArc is recommended to be between +/- 5 deg \n";
				countArcDynCollAngleRemarks++;
			}
			return cResults;
		}


	}
}
