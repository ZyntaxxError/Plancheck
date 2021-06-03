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
using System.Collections;
using System.Text;
using System.Diagnostics;



/* TODO: 
 * TODO: Clean up extensions, some methods are not actual extensions
 * TODO Plan sum; foreach plan: check, easier to do in build and wpf...
 * TODO: For dynamic plans; check if verification plan exists only if plan status is planning approved 
 * TODO: z-value for imaging, recommend extended CBCT
 * TODO: (double[2]) var gRange = plan.GetGantryRange() to check for possible collission
 * collision at beam or traveling to beam, assuming starts at gantry 0? imaging?
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

		public enum UnitLength
		{
			mm = 1,
			cm = 10,
			dm = 100,
			m = 1000
		}

		public enum Direction
		{
			Cranial,
			Caudal,
			None
		}







		public static double ScaleConv(double val, Scale fromScale, Scale toScale, Machine.Axis axis, bool ExtendedFlag = false)
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
									if (val >= 0 && val <= 180 && ExtendedFlag == false)
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



	public static class Machine
	{
		public const int SID = 1000;  // Source-isocenter-distance in mm
									  // speeds and accelerations (except Collimator) are mean values from trajectory logs, parameters vary somewhat depending on acceleration or deceleration, gravity etc 
									  // but should be close enough to calculate a fair estimation of beam-on time. Collimator speed taken from TrueBeam Specs, collimator acceleration is arbitrarily chosen.
		public const int GantryMaxSpeed = 6;    // deg/s
		public const int CollMaxSpeed = 15;     // deg/s,  2.5 RPM with no accessory (1 RPM with accessory)
		public const float JawYMaxSpeed = 22.5f;
		public const float JawXMaxSpeed = 22.5f;
		public const int GantryMaxAcc = 16;     // deg/s^2 
		public const int CollMaxAcc = 60;       // TODO: TBD
		public const int JawYMaxAcc = 50;       // mm/s^2
		public const int JawXMaxAcc = 160;      // mm/s^2

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
				{ "6X-FFF", 1400 },
				{ "10X-FFF", 2400 },
				{ "6E", 1000 },
				{ "9E", 1000 },
				{ "12E", 1000 }
			};
	}

	public static class PlanExtensions
	{
		public static bool AutomationPrerequisites(this PlanSetup plan)
		{
			//TODO: Missing check for same accessories for all beams. Define accessories... wedges included?
			if (plan.SingleIso() && plan.SingleCouchPosition() && plan.SingleEnergy())
			{
				return true;
			}
			else
			{
				return false;
			}
		}

		public static bool SingleIso(this PlanSetup plan)
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

		public static bool SingleEnergy(this PlanSetup plan)
		{
			bool sEnergy = true;
			Beam firstBeam = plan.Beams.Where(b => !b.IsSetupField).OrderBy(b => b.Id).First();
			foreach (var beam in plan.Beams.Where(b => !b.IsSetupField).OrderBy(b => b.Id))
            {
                if (beam.EnergyModeDisplayName.Equals(firstBeam.EnergyModeDisplayName))
                {
                    continue;
                }
                else
                {
					sEnergy = false;
					break;
				}
            }
            return sEnergy;
		}

		public static bool SingleCouchPosition(this PlanSetup plan)
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


		public static bool AnyWedges(this PlanSetup plan)
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
	}

	public static class BeamExtensions
	{

		public static bool MaxDoseRateUsed(this Beam beam)
		{
			int maxDoseRate = 0;
			if (!Machine.DoseRates.TryGetValue(beam.EnergyModeDisplayName, out maxDoseRate))
			{
				MessageBox.Show($"{beam.EnergyModeDisplayName} is not recognised as an active energy mode!");
				return true; // Ugly code to prevent it from crashing, fix this
			}

			if (beam.DoseRate == maxDoseRate)
			{
				return true;
			}
			else
			{
				return false;
			}
		}


		/// <summary>
		/// Calculates an estimated beam on time from control points in beam, works for dynamic as well as static fields. TODO: EDW unhandled, underestimates 
		/// time for wedges, especially large fields (Y) with small number of MU, need to get the STT-tables and calculate this separately 
		/// </summary>
		/// <param name="beam"></param>
		/// <returns></returns>
		public static int EstimatedBeamOnTime(this Beam beam)
		{
			double timeOffset = 0.9; // s, added time for startup beam stabilisation 0.4 s, empirical estimation from trajectory logs. Stopping time last control point aproximated to 0.5 s
			double time = 0;
			int maxDoseRate = beam.DoseRate / 60; // MU/s
			double totalMU = beam.Meterset.Value;
			double deltaMU;

			int nrOfJaws = 4;

			int nrCP = beam.ControlPoints.Count;
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

			for (int i = 1; i < beam.ControlPoints.Count; i++)
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
				// this might be completely off... problably need to study how the control system acts between control points, but closer
				// to the truth than shown in control points in Eclipse

				gantrySpeed[i] = deltaGantry[i] / cpTime[i];

				for (int j = 0; j < nrOfJaws; j++)
				{
					jawSpeed[i, j] = deltaJaw[i, j] / cpTime[i]; // mm/s
				}

				time += cpTime[i];

				//debug list of control points for comparison with trajectory logs.
				controlPointList.AppendLine(i + "\t" + cpTime[i].ToString("0.00") + "\t" + "\tGspeed:" + gantrySpeed[i].ToString("0.0"));
				for (int j = 0; j < nrOfJaws; j++)
				{
					controlPointList.Append((jawSpeed[i, j] / 10).ToString("0.0") + "\t"); // Convert to cm/s to compare with log files
				}
				doseRate[i] = deltaMU * 60 / cpTime[i];     // calculation of doserate to compare with log files MU/min
				controlPointList.Append("\t" + doseRate[i]);
			}

			string cpInfo = controlPointList.ToString();

			//MessageBox.Show(cpInfo);

			return (int)(time + timeOffset);
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

			//The sign of s and vmax (acceleration?) should always be the same when calculating minimum time 
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
		/// TODO: Extended range for gantry angles unhandled, no flag for this in esapi (actually there is in 16.1.. wohoo!). 
		/// </summary>
		/// <param name="angleStart"></param>
		/// <param name="angleEnd"></param>
		/// <returns></returns>
		public static double DeltaAngle(Double angleStart, Double angleEnd, Conversion.Scale scale, Machine.Axis axis)
		{
			//Convert to Varian machine scale before calculating delta angle 
			angleStart = Conversion.ScaleConv(angleStart, scale, Conversion.Scale.VarianIEC, axis);
			angleEnd = Conversion.ScaleConv(angleEnd, scale, Conversion.Scale.VarianIEC, axis);
			return Math.Abs(angleEnd - angleStart);
		}
	}


	public class StructureModel
    {
		/*
 * Order: by Id:  PTV, E_PTV, ITV, CTV, GTV, BODY, BolusIncluded in calc, Support,	 avoidance, organs and the rest, help structures, Bolus not included,	helpBolus, empty structures
 * 
 * Dicom type	description		volume		nrofparts	assigned hu		comment in eclipse	
 */
		public string Id { get; set; }
        public string Type { get; set; }
        public string Description { get; set; }
        public bool Approved { get; set; }
        public int? DoseValueInId { get; set; }
        public double VolumeCC { get; set; }
        public int NrOfParts { get; set; }
        public double? AssignedHU { get; set; }
        public string AriaComment { get; set; }
        public List<string> Warnings { get; set; }
		public List<string> Errors { get; set; }
		public List<string> Comments { get; set; }
        public bool HighRes { get; set; }

        public enum  StructureType
        {
			PTV,
			ITV,
			CTV,
			GTV,
			EvalPTV,
			BODY,
			Support,

        }

		//Dicom type: Possible values are "AVOIDANCE", "CAVITY", "CONTRAST_AGENT", "CTV", "EXTERNAL", "GTV", "IRRAD_VOLUME", 
		//"ORGAN", "PTV", "TREATED_VOLUME", "SUPPORT", "FIXATION", "CONTROL", and "DOSE_REGION". 
		// where is ITV? 


	}

	public class StructureViewModel
    {

		public StructureViewModel(StructureSet ss)
        {
        //List<StructureModel> Structures = new List<StructureModel>();
			GetStructureData(ss);
		}

		public List<StructureModel> Structures = new List<StructureModel>();

		// structures included resp excluded in calculation (outside body and not bolus, empty)
		// sort order type; ptv, ctv, itv, gtv 

		private void GetStructureData(StructureSet ss)
        {

			List<Structure> ActiveStructures = new List<Structure>();
			List<Structure> EmptyStructures = new List<Structure>();
			string debug = string.Empty;
			//EmptyStructures = ss.Structures.Where(s => s.HasSegment == false).ToList();
			try
			{
				ActiveStructures = ss.Structures.Where(s => s.HasSegment).OrderBy(s => s.Id).ToList();
				foreach (Structure str in ActiveStructures)
				{
					Structures.Add(new StructureModel
					{
						Warnings = new List<string>(),
						Errors = new List<string>(),
						Comments = new List<string>(),
						Id = str.Id,
						Type = str.DicomType,
						Description = str.StructureCode.Code.FirstOrDefault().ToString(),
						Approved = str.IsApproved,
						AssignedHU = GetAssignedHU(str),
						VolumeCC = str.Volume,
						NrOfParts = str.GetNumberOfSeparateParts()




				});


				}



				// add comments, warnings and errors
				foreach (var strModel in Structures)
				{
					if (strModel.AssignedHU.HasValue && strModel.Id.Contains("Couch") == false)
					{
						strModel.Comments.Add("Houndsfieldvärde overrajdat");
					}
					debug += strModel.Id + "OK\n";




















				}
			}
            catch (Exception e)
            {
				MessageBox.Show(e.Message + debug);
            }
        }


		public double? GetAssignedHU(Structure str)
        {
			double? hu = null;
			double huu;
            if (str.GetAssignedHU(out huu))
            {
				hu = huu;
            }
			return hu;
        }







	}









	public static class StructureExtensions
	{
		/// <summary>
		/// Calculates geometric center of structure bounds in Dicom coordinates
		/// TODO: Check if this is dependent on patient orientation, probably not, but check it anyway
		/// </summary>
		/// <param name="structure"></param>
		/// <returns></returns>
		public static VVector GeometricCenter(this Structure structure)
		{
			VVector geomCenter = new VVector
			{
				x = structure.MeshGeometry.Bounds.X + structure.MeshGeometry.Bounds.SizeX / 2,
				y = structure.MeshGeometry.Bounds.Y + structure.MeshGeometry.Bounds.SizeY / 2,
				z = structure.MeshGeometry.Bounds.Z + structure.MeshGeometry.Bounds.SizeZ / 2
			};
			return geomCenter;
		}
		/// <summary>
		/// Calculates geometric center of structure in the closest slice to given position
		/// </summary>
		/// <param name="structure"></param>
		/// <returns></returns>
		public static VVector GeometricCenterClosestSlice(this Structure structure, StructureSet ss, VVector pos)
		{
			Image image = ss.Image;
			int slice = image.ClosestSlice(pos);
			VVector[][] contour = structure.GetContoursOnImagePlane(slice);
			// remake the jagged array to 1-D List
			List<VVector> flattenPoints = new List<VVector>();
			foreach (var arr in contour)
			{
				foreach (var point in arr)
				{
					flattenPoints.Add(point);
				}
			}

			VVector geomCenter = new VVector
			{
				x = (flattenPoints.OrderBy(p => p.x).First().x + flattenPoints.OrderBy(p => p.x).Last().x) / 2,
				y = (flattenPoints.OrderBy(p => p.y).First().y + flattenPoints.OrderBy(p => p.y).Last().y) / 2,
				z = (flattenPoints.OrderBy(p => p.z).First().z + flattenPoints.OrderBy(p => p.z).Last().z) / 2
			};
			return geomCenter;
		}




		/// <summary>
		/// Maximum distance [mm] in x-y-plane from given Dicom coordinate to exitpoint of structure
		/// </summary>
		/// <param name="structure"></param>
		/// <param name="isoPos"></param>
		/// <returns>double; max distance in mm</returns>
		public static double MaxDistToStructureExit(this Structure structure, VVector pos)
		{
			// rotate around the structure with origo in iso and find the maximum distance
			double maxDistance = 0;
			// use outer bounds to determine necessary search radius
			double searchRadius = structure.MaxDistToBoundsExit(pos);
			int samples = (int)searchRadius;
			for (int angle = 0; angle < 360; angle++)
			{
				double x = pos.x + searchRadius * Math.Cos(angle * Math.PI / 180);
				double y = pos.y + searchRadius * Math.Sin(angle * Math.PI / 180);
				VVector endPoint = pos;
				endPoint.x = x;
				endPoint.y = y;

				// Only check the distance if the profile actually intersects with the structure
				var testIntersect = structure.GetSegmentProfile(pos, endPoint, new BitArray(samples)).Where(s => s.Value == true);
				if (testIntersect.Any())
				{
					SegmentProfilePoint structureExitPoint = testIntersect.Last();

					//convert from Type Point to VVector, may be unneccesary unless position is of interest
					VVector exitPos = new VVector
					{
						x = structureExitPoint.Position.x,
						y = structureExitPoint.Position.y,
						z = structureExitPoint.Position.z
					};

					double distance = Math.Sqrt(Math.Pow(pos.x - exitPos.x, 2) + Math.Pow(pos.y - exitPos.y, 2));
					if (distance > maxDistance)
					{
						maxDistance = distance;
					}
				}
			}
			return maxDistance;
		}

















		/// <summary>
		/// Maximum distance [mm] in x-y-plane from given Dicom coordinate to exitpoint of structure
		/// TODO: investigate if this could be used instead: VVector[][] contour = structure.GetContoursOnImagePlane(slize_z)
		/// </summary>
		/// <param name="structure"></param>
		/// <param name="isoPos"></param>
		/// <returns>double; max distance in mm</returns>
		public static double MaxDistToStructureContourExit(this Structure structure, StructureSet ss, VVector pos)
		{
			Image image = ss.Image;
			int slice = image.ClosestSlice(pos);
			double maxDistance = 0;
			VVector[][] contour = structure.GetContoursOnImagePlane(slice);

			foreach (var arr in contour)
			{
				foreach (var point in arr)
				{
					double distance = Math.Sqrt(Math.Pow(pos.x - point.x, 2) + Math.Pow(pos.y - point.y, 2));

					if (distance > maxDistance)
					{
						maxDistance = distance;
					}
				}
			}
			return maxDistance;
		}









		/// <summary>
		/// Maximum distance [mm] in x-y-plane from given Dicom coordinate to outer bounds of structure.
		/// Primarily used to get a search radius for structure in a x-y-plane
		/// </summary>
		/// <param name="structure"></param>
		/// <param name="pos"></param>
		/// <returns></returns>
		public static double MaxDistToBoundsExit(this Structure structure, VVector pos)
		{
			double lat1 = structure.MeshGeometry.Bounds.X;
			double lat2 = structure.MeshGeometry.Bounds.X + structure.MeshGeometry.Bounds.SizeX;
			/*
			if (image.ImagingOrientation == PatientOrientation.HeadFirstSupine)
            {
				lat2 = structure.MeshGeometry.Bounds.X + structure.MeshGeometry.Bounds.SizeX;
            }
            else
            {
				lat2 = structure.MeshGeometry.Bounds.X - structure.MeshGeometry.Bounds.SizeX;
			}
			*/
			// maybe use target.MeshGeometry.Positions.Max(p => p.Z)
			double vrt1 = structure.MeshGeometry.Bounds.Y;
			double vrt2 = structure.MeshGeometry.Bounds.Y + structure.MeshGeometry.Bounds.SizeY;

			double[] dist = new double[4];

			dist[0] = Math.Sqrt(Math.Pow(pos.x - lat1, 2) + Math.Pow(pos.y - vrt2, 2));
			dist[1] = Math.Sqrt(Math.Pow(pos.x - lat2, 2) + Math.Pow(pos.y - vrt2, 2));
			dist[2] = Math.Sqrt(Math.Pow(pos.x - lat1, 2) + Math.Pow(pos.y - vrt1, 2));
			dist[3] = Math.Sqrt(Math.Pow(pos.x - lat2, 2) + Math.Pow(pos.y - vrt1, 2));

			return dist.Max();
		}


	}


	public static class ImageExtensions
	{
		/// <summary>
		/// Calculates closest slice number in image from given Dicom position.
		/// </summary>
		/// <param name="image"></param>
		/// <param name="pos"></param>
		/// <returns>Slice number enumerated in caudal-cranial direction</returns>
		public static int ClosestSlice(this Image image, VVector pos)
		{
			double distanceFromOrigin = Math.Abs(pos.z - image.Origin.z) / image.ZRes;
			int sliceNumber = (int)Math.Round(distanceFromOrigin) + 1;
			return sliceNumber;
		}
		public static int ClosestSlice(this Image image, double zPos)
		{
			double distanceFromOrigin = Math.Abs(zPos - image.Origin.z) / image.ZRes;
			int sliceNumber = (int)Math.Round(distanceFromOrigin) + 1;
			return sliceNumber;
		}

	}




	//******************************************************************   SCRIPT START *************************************************


	// TODO: bool IsGantryExtended
	// TODO: change "structure with assigned HU" when checking calculated plan to; "structures with assigned hu included in calculation"
	//: TODO, make functions for TMI OPT work without plan, should be possible to recommend isocenters (even cm) and limits of last z_ptv
	// based on ID of image (series properties?) Prereq: BODY, origo and couch ok
	//DEBUG
	//using System.Diagnostics;
	//var sw = Stopwatch.StartNew();
	//Console.WriteLine("Tid " + sw.ElapsedMilliseconds + " ms:  ");


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
			TBI,
			TMI
		}

		public void Execute(ScriptContext context)
		{
			if ((context.PlanSum == null && context.PlanSetup == null) || context.StructureSet == null)
			{
				if (context.StructureSet == null)
				{
					MessageBox.Show("Please select a plan in the active context window.");
				}
				else
				{
					StructureSet ss = context.StructureSet;


					var structView = new StructureViewModel(ss);

					string structMessage = string.Empty;

                    foreach (var item in structView.Structures)
                    {
						structMessage += item.Id + "\t\t" + item.Type + "\n";
                        if (item.Comments.Count > 0)
                        {
							structMessage += item.Comments[0];
						}
                        if (item.AssignedHU.HasValue)
                        {
							structMessage += item.AssignedHU.Value.ToString("0.0");
						}
					}
					MessageBox.Show(structMessage);

					


					string messageTitle = "Quick check on structure set" + ss.Id;
					string message = CheckStructureSet(ss) +
					CheckCouchStructure(ss);
					MessageBox.Show(message, messageTitle);
				}
			}
			else
			{
				if (context.PlanSum != null)
				{
					PlanSum psum = context.PlanSum;
					CheckPlanSum(psum);  //TODO: check plan category of included plans in plan sum
				}
				else
				{
					string usr = context.CurrentUser.Id;
					PlanSetup plan = context.PlanSetup;
					//string planCreator = plan.HistoryUserName.ToString();
					PlanCat planCat = PlanCat.Unknown;
					CheckPlanCategory(plan, ref planCat);
					List<string> relevantDocuments = new List<string>();
					Course course = context.Course;
					string courseIntent = context.Course.Intent;
					string courseId = context.Course.Id;
					string messageTitle = "Quick check on " + courseId + " " + plan.Id;
					string message = string.Empty;

					message +=
					"Hello " + usr + "\n" +
					"Assumed plan category: " + planCat + "\n\n" +
					CheckCourseIntent(courseIntent, plan, planCat) + "\n" +
					CheckPlanNamingConvention(plan) +
					CheckClinProt(plan) +
					CheckTargetVolumeID(plan, plan.StructureSet) +
					SBRTCtvPtvPlanNumbering(plan, planCat) +
					CheckForBolus(plan, ref relevantDocuments) +
					CheckFieldNamingConvention(course, plan);


					if (planCat == PlanCat.Electron)
					{
						message += CheckElectronPlan(plan);
					}
					else
					{
						switch (planCat)
						{
							case PlanCat.Unknown:
								break;
							case PlanCat.Electron:
								break;
							case PlanCat.SBRT:

								break;
							case PlanCat.TBI:
								message += CheckBodyCenter(plan);
								break;
							case PlanCat.TMI:
								if (plan.Id.Contains("OPT"))
								{
									message += CheckBodyCenter(plan);
									message += TMIGetMaxSSDJunctions(plan);
								}
								break;
							default:
								break;
						}

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

					MessageBox.Show(message, messageTitle);

				}
			}
		}



		#region TMI operations

		//TODO: compare body bounds with image size and check if margin is enough
		// assuming we have two CT, BODY is contoured on both, possible to calculate total length; detemining factor for iso distances and field sizes
		// Dicom coordinates last iso HFS to FFS-CT?


		//TODO: error handling; check for prereq: Structures needed etc
		// for FFS AND Base dose plan: approved registration with only translations (registration direction?)
		// check naming of image and structure set




		// Optimizer checks for overlapping structures with support structures; ex: bolus and couch




		//     Structure structure = ss.Structures.Where(s => s.Id == "PTV_Total").SingleOrDefault();
		//     int margin = 0; // 10 mm margin added in critical places if no PTV_Total contoured
		//if (structure == null || structure.IsEmpty)
		//{
		//	Structure structure = ss.Structures.Where(s => s.Id == "BODY").SingleOrDefault();
		//     margin = 10;
		// }


		/// <summary>
		/// Gets the maximum SSD for each junction in a TBI/TMI plan
		/// ASSUMPTION: same fixed field size for every iso giving the junctions in the middle between the isos
		/// </summary>
		/// <param name="plan"></param>
		/// <returns></returns>
		private string TMIGetMaxSSDJunctions(PlanSetup plan)
		{
			// need to check if all isocenters are positioned in the same lat and vrt value
			StructureSet ss = plan.StructureSet;
			Image image = ss.Image;
			Structure structure = ss.Structures.Where(s => s.Id == "PTV_Total").SingleOrDefault();  // TODO: better to take dicom type, might want to use PTV_Total instead
			string junkPosZ = string.Empty;
			// get list of isocenter positions
			List<VVector> isos = TMIIsocentersInPlan(plan);





			// debug
			for (int i = 0; i < isos.Count; i++)
			{
				VVector junkEclipse = image.DicomToUser(isos[i], plan);
				junkPosZ += "iso " + i + ":\t" + (junkEclipse.z / 10).ToString("0.0") + "\n";
				//public IEnumerable<Registration> Registrations { get; }

			}

			if (isos.Count == 1)
			{
				junkPosZ += "\n** Only one isocenter, no checks can be done for this\n";
			}
			else if (isos.Any(i => i.x != isos[0].x || i.y != isos[0].y))
			{
				junkPosZ += "\n** Please adjust all isocenters to identical positions in vrt and lat\n";
			}

			List<VVector> junctionPositions = new List<VVector>(); // TODO: assumes all field sizes the same
			VVector junction = new VVector();
			for (int i = 1; i < isos.Count; i++)
			{
				junction.x = isos[i].x;
				junction.y = isos[i].y;
				junction.z = (isos[i].z + isos[i - 1].z) / 2;
				junctionPositions.Add(junction);
				VVector junkEclipse = image.DicomToUser(junction, plan);
				junkPosZ += "junktion " + i + ":\t" + (junkEclipse.z / 10).ToString("0.0") + "\n";
			}

			double searchRadius = structure.MaxDistToBoundsExit(junction);
			VVector exitPos = new VVector();
			VVector[] maxDistancePosition = new VVector[junctionPositions.Count];
			double[] maxDistance = new double[junctionPositions.Count];
			double[] maxDistance2 = new double[junctionPositions.Count];
			double[] maxDistanceUp = new double[junctionPositions.Count];
			double[] maxDistanceDown = new double[junctionPositions.Count];
			maxDistance[0] = 0.0;
			maxDistance2[0] = 0.0;
			maxDistanceUp[0] = 0.0;
			maxDistanceDown[0] = 0.0;

			for (int j = 0; j < junctionPositions.Count; j++)
			{
				maxDistance2[j] = structure.MaxDistToStructureContourExit(ss, junctionPositions[j]);
				VVector Up = junctionPositions[j];
				VVector Down = junctionPositions[j];
				Up.z += image.ZRes;
				Down.z -= image.ZRes;
				maxDistanceUp[j] = structure.MaxDistToStructureContourExit(ss, Up);
				maxDistanceDown[j] = structure.MaxDistToStructureContourExit(ss, Down);






				// rotate around the structure with origo in junction and find the maximum distance
				for (int angle = 0; angle < 360; angle++)
				{
					double x = junction.x + searchRadius * Math.Cos(angle * Math.PI / 180);
					double y = junction.y + searchRadius * Math.Sin(angle * Math.PI / 180);
					VVector endPoint = junctionPositions[j];
					endPoint.x = x;
					endPoint.y = y;


					// condition if segment profile don't intersect with the structure. 
					var testIntersect = structure.GetSegmentProfile(junctionPositions[j], endPoint, new BitArray(200)).Where(s => s.Value == true);
					if (testIntersect.Count() > 0)
					{


						SegmentProfilePoint structureExitPoint = structure.GetSegmentProfile(junctionPositions[j], endPoint, new BitArray(200)).Where(s => s.Value == true).Last();

						//convert from Type Point to VVector
						exitPos.x = structureExitPoint.Position.x;
						exitPos.y = structureExitPoint.Position.y;
						exitPos.z = structureExitPoint.Position.z;

						double distance = Math.Sqrt(Math.Pow(junctionPositions[j].x - exitPos.x, 2) + Math.Pow(junctionPositions[j].y - exitPos.y, 2));

						if (distance > maxDistance[j])
						{
							maxDistance[j] = distance;
							maxDistancePosition[j] = exitPos;
						}
					}

				}
			}
			// Note that all positions are still in dicom coordinates and distances in mm

			double fz = 0;
			if (plan.TreatmentOrientation == PatientOrientation.HeadFirstSupine)
			{
				fz = 13.5;
			}
			else
			{
				fz = 14.0;
			}

			for (int i = 0; i < junctionPositions.Count; i++)
			{
				VVector maxDistPosEclipse = image.DicomToUser(maxDistancePosition[i], plan);
				junkPosZ += "\n\nmaxpos(x,y,z) " + ":\t(" + (maxDistPosEclipse.x / 10).ToString("0.0") + ", " + (maxDistPosEclipse.y / 10).ToString("0.0") + ", " + (maxDistPosEclipse.z / 10).ToString("0.0") + ")\n";
				junkPosZ += "Max distance:\t" + (maxDistance[i] / 10).ToString("0.0") + "\n";
				junkPosZ += "Max distance2:\t" + (maxDistance2[i] / 10).ToString("0.0") + "\n";
				junkPosZ += "Max distanceUp:\t" + (maxDistanceUp[i] / 10).ToString("0.0") + "\n";
				junkPosZ += "Max distanceDown:\t" + (maxDistanceDown[i] / 10).ToString("0.0") + "\n";
				junkPosZ += "Max delta to cover the junction in all angles: \t" + (2 * fz * (100 - (maxDistance[i] / 10)) / 100).ToString("0.0") + "\n";
			}

			junkPosZ += "\nMin SSD:\t" + (100 - maxDistance.Max() / 10).ToString("0.0") + "\n\n";

			if (plan.TreatmentOrientation == PatientOrientation.HeadFirstSupine)
			{
				junkPosZ += TMICranialCoverage(plan);
			}



			List<VVector> tmiRecom = new List<VVector>();
			if (plan.TreatmentOrientation == PatientOrientation.HeadFirstSupine)
			{
				tmiRecom = TMISuggestedIsocentersHFS(plan);

			}
			else
			{
				tmiRecom = TMISuggestedIsocentersFFS(plan);
			}


			TMIcheckIfCropped(plan);


			List<VVector> tmiRecomEclipse = new List<VVector>();
			foreach (var iso in tmiRecom)
			{
				tmiRecomEclipse.Add(image.DicomToUser(iso, plan));
			}




			string recom = "\n\n Suggested isocenters: (x,y,z) \n";

			foreach (var i in tmiRecomEclipse)
			{
				recom += (i.x / 10).ToString("0.0") + "\t" + (i.y / 10).ToString("0.0") + "\t" + (i.z / 10).ToString("0.0") + " \n";
			}

			// Recommend position for Z_PTV5

			if (plan.TreatmentOrientation == PatientOrientation.HeadFirstSupine)
			{
				recom += "\nRecommended position for Z_PTV5 is between isocenter 5 and 6.";
			}

			MessageBox.Show(recom);

			return junkPosZ;


			//Delta = 2*x*SSD/iso to overlap, where x is jaw position 


		}

		/*

)[EditorBrowsableAttribute(EditorBrowsableState.Advanced)]
public PlanSetup CopyPlanSetup(
	PlanSetup sourcePlan,
	Image targetImage,
	Registration registration,
	StringBuilder outputDiagnostics
)


		 */
		// TODO check position of couch and check how algorithms handles overlap of couch with body, PTV and bolus; OK optimizer checks this
		// marker on HFS "Z_92" 92 cm below cranial part of necksupport intendention
		// register image with CX DP FFS, transfer structure to registered image, approve registration and marker
		// script checks minimum delta needed for HFS, calculates position of last iso HFS in "CX DP FFS" using marker as reference 
		// When calculating FFS check if it is OK to separate delta, and in that case allow smaller deltashift for HFS -> more overlap means prob better plan
		// deltashift in steps of 5 mm would make planning easier...
		// Order for manual plan;  suggested iso, place isos and copy to FFS i.e make base dose plan in HFS using image in FFS, check suggested isos, if mucho overlap, decrease delta HFS

		// operlap structures in border between HFS and FFS to 


		// if minimum overlap > 4 * 1 cm && delta == 20 -> decrease delta for hfs (check that we don't add isocenters in FFS)
		// public List<VVector> TMISuggestedIsocentersHFS(PlanSetup plan, int maximumDelta = 200;)
		public List<VVector> TMISuggestedIsocentersHFS(PlanSetup plan)
		{
			StructureSet ss = plan.StructureSet;
			Image image = plan.StructureSet.Image;
			List<VVector> isos = new List<VVector>();
			double fieldSize = 135.0;
			int maximumDelta = 200;
			int minimumDelta = 180; // minimum acceptable delta, should probably try to avoid this however
			int deltaStep = 10; // 10 mm step between investigated delta shifts
			int cranialMargin = 5; // Extra cranial margin for buildup, should be enough as the isocenter is determined by the z-position which covers all angles
			int latShiftLimit = 5; // limit [mm] for when to shift the isocenters from origo position in lateral direction
			int minOverlap = 0;  // mm overlap where minimum occurs, i.e. maximum distance from isocenter plane 


			Structure ptvTotal = ss.Structures.Where(s => s.Id == "PTV_Total").SingleOrDefault();  // TODO: better to take dicom type
			VVector ptvGeoCenter = ptvTotal.GeometricCenter();

			Structure lungTotal = ss.Structures.Where(s => s.Id == "LungTotal" || s.Id == "Lung Total").SingleOrDefault();  // TODO: better to take dicom type
			VVector lungGeoCenter = lungTotal.GeometricCenter();

			// Set temporary position for first cranial isocenter as user origin
			VVector firstIsoCran = ss.Image.UserOrigin;

			// Vrt resonable position is estimated by weighting geometric center of lungs and PTV_Total
			firstIsoCran.y = (lungGeoCenter.y * 60 + ptvGeoCenter.y * 40) / 100;

			// Long position for first isocenter is determined by the z value where the given field size covers ptv in all angles
			firstIsoCran = ZposCoverAllAngles(ptvTotal, firstIsoCran, fieldSize, true);
			firstIsoCran.z += cranialMargin;

			// Round vrt to nearest cm and long to nearest mm from user origo (lat is per def zero) 
			firstIsoCran.y = RoundToUser(plan, firstIsoCran, Conversion.UnitLength.cm).y;
			firstIsoCran.z = RoundToUser(plan, firstIsoCran).z;

			List<VVector> junctionPositions = new List<VVector>();

			// Iteration; check if overlap occurs in the junction regions and their vicinity (+/- 5 mm) for all angles, if not decrease delta
			// and check again until minimum acceptable delta reached. 
			for (int delta = maximumDelta; delta >= minimumDelta; delta -= deltaStep)
			{
				isos.Add(firstIsoCran);
				VVector iso = firstIsoCran;
				for (int i = 1; i < 6; i++)
				{
					iso.z -= delta;
					isos.Add(iso);
				}

				// preliminary isocenters based on first isocenter and delta shift, make list of junctionpositions to check lateral shift and overlap
				VVector junction = new VVector();
				for (int i = 1; i < isos.Count; i++)
				{
					junction = isos[i];
					junction.z = (isos[i].z + isos[i - 1].z) / 2;
					junctionPositions.Add(junction);
				}

				// check if ptv in junctions (closest slice) and lungsTotal is displaced relative user origin in the lateral direction
				// Ignore first junction as this might be located in neck/shoulder area where lat might shift rapidly, replace it with first iso.
				// If shift is larger than the set limit; set new lat pos.
				// First iso long should still be ok even with lat shift of a couple of cm.
				double optLat = lungGeoCenter.x + ptvTotal.GeometricCenterClosestSlice(ss, isos[0]).x;
				for (int i = 1; i < junctionPositions.Count; i++)
				{
					optLat += ptvTotal.GeometricCenterClosestSlice(ss, junctionPositions[i]).x;
				}
				optLat /= (junctionPositions.Count + 2);
				optLat = (int)Math.Round(optLat, 0);

				if (Math.Abs(optLat - image.UserOrigin.x) > latShiftLimit)
				{
					for (int i = 0; i < isos.Count; i++)
					{
						VVector temp = isos[i];
						temp.x = optLat;
						isos[i] = temp;
					}

					for (int i = 0; i < junctionPositions.Count; i++)
					{
						VVector temp = junctionPositions[i];
						temp.x = optLat;
						junctionPositions[i] = temp;
					}
				}

				//find maximum distance from each junction position (lng) in isocenter plane (lat, vrt) to structure exitpoint (i.e. minimum SSD for full arc)
				double maxDistance = 0;
				for (int i = 0; i < junctionPositions.Count; i++)
				{
					VVector pos = junctionPositions[i];
					double dist = ptvTotal.MaxDistToStructureExit(pos);
					if (dist > maxDistance)
					{
						maxDistance = dist;
					}
				}

				// break if condition is satisfied, otherwise reduce delta and check again
				if (delta + minOverlap < 2 * fieldSize * (Machine.SID - (maxDistance)) / Machine.SID || delta == minimumDelta)
				{
					break;
				}
				else
				{
					isos.Clear();
					junctionPositions.Clear();
				}
			}
			return isos;
		}






		/// <summary>
		/// Rounds a dicom position to nearest even unit in the user coordinate system (i.e relative to user origin)
		/// </summary>
		/// <param name="plan"></param>
		/// <param name="pos"></param>
		/// <param name="unit"></param>
		/// <returns>position in dicom coordinates</returns>
		public VVector RoundToUser(PlanSetup plan, VVector pos, Conversion.UnitLength unit = Conversion.UnitLength.mm)
		{

			int unitConv = (int)unit;
			Image image = plan.StructureSet.Image;
			VVector eclipseRound = image.DicomToUser(pos, plan);
			eclipseRound.x = (int)Math.Round(eclipseRound.x / unitConv, 0) * unitConv;
			eclipseRound.y = (int)Math.Round(eclipseRound.y / unitConv, 0) * unitConv;
			eclipseRound.z = (int)Math.Round(eclipseRound.z / unitConv, 0) * unitConv;

			return image.UserToDicom(eclipseRound, plan);
		}

		// TODO: what happens if enough with image in HFS
		// prereq: user origo placed in iso 5 from HFS-plan
		// can use coordinates from Base dose HFS, as the isocenters may not have been placed according to suggestion... 
		// TODO: Special case that need to be handled; short patient, enough with 1 isocenters, well... no need for script in that case if not used for auto plan
		public List<VVector> TMISuggestedIsocentersFFS(PlanSetup plan)
		{
			StructureSet ss = plan.StructureSet;
			List<VVector> isos = new List<VVector>();
			double fieldSizeFFS = 140.0;
			double fieldSizeHFS = 135.0;
			int maxDeltaHFS = 200;
			int maxDeltaFFS = 240;
			int minDeltaFFS = 180;
			int deltaStep = 10; // 10 mm step between investigated delta shifts
			Image image = ss.Image;
			int caudalMargin = 1; // Extra caudal margin for buildup, this is minimum but almost always increases as calulated delta is rounded up to nearest cm
			int minOverlap = 0;  // mm overlap where minimum occurs, i.e. maximum distance from isocenter plane 

			Structure ptvTotal = ss.Structures.Where(s => s.Id == "PTV_Total").SingleOrDefault();  // TODO: better to take dicom type

			// origo should be placed at the last isocenter in HFS plan (P5)
			VVector lastIsoHFS = ss.Image.UserOrigin;

			// Determine position of first iso FFS separately as field size for HFS is smaller which means junction in different place
			// Delta coordinates will be different and will be sent from origo, want reassuring overlap so calculate with field size from HFS
			VVector firstIsoFFS = lastIsoHFS;
			for (int delta = maxDeltaHFS; delta >= minDeltaFFS; delta -= deltaStep)
			{
				firstIsoFFS.z = lastIsoHFS.z - delta;
				VVector junction = lastIsoHFS;
				junction.z = (firstIsoFFS.z + lastIsoHFS.z) / 2;
				double maxDistance = ptvTotal.MaxDistToStructureExit(junction);

				if (delta + minOverlap < 2 * fieldSizeHFS * (Machine.SID - (maxDistance)) / Machine.SID || delta == minDeltaFFS)
				{
					break;
				}
			}

			// Preliminary position for last iso FFS is determined by scanning for the position where the given field size covers the target in all angles
			VVector lastIsoFFS = ZposCoverAllAngles(ptvTotal, firstIsoFFS, fieldSizeFFS, false);
			lastIsoFFS.z -= caudalMargin;

			// start value for number of isos, i.e. minimum number of isocenters using maximum delta
			int minNrOfIsos = (int)Math.Ceiling(Math.Abs(firstIsoFFS.z - lastIsoFFS.z) / maxDeltaFFS);

			List<VVector> junctionPositions = new List<VVector>();

			// Iteration; start with maximum delta, calculate nr of isos necessary, check if overlap occurs in the 
			//junction regions for all angles, if not; increase nr of isos and check again until minimmum acceptable delta reached. 
			for (int nriso = minNrOfIsos; nriso <= 8; nriso++)
			{
				isos.Add(firstIsoFFS);
				VVector iso = firstIsoFFS;


				// change delta to even cm based on number of isos and back to mm for calculation
				double firstToLastIsoMm = Math.Abs(firstIsoFFS.z - lastIsoFFS.z); // mm
				int delta = (int)Math.Ceiling(firstToLastIsoMm / (10 * nriso)) * 10;

				//MessageBox.Show("FirstToLastIsoMM: " + firstToLastIsoMm + "\tDelta  :"  + delta + "\tNr of iso: " + nriso + "\n");

				for (int i = 1; i <= nriso; i++)
				{
					iso.z -= delta;
					isos.Add(iso);
				}

				// preliminary isocenters based on first isocenter and delta shift, make list of junctionpositions to check overlap
				VVector junction = new VVector();
				for (int i = 1; i < isos.Count; i++)
				{
					junction = isos[i];
					junction.z = (isos[i].z + isos[i - 1].z) / 2;
					junctionPositions.Add(junction);
				}

				//find maximum distance from each junction position (lng) in isocenter plane (lat, vrt) to structure exitpoint (i.e. minimum SSD)
				double maxDistance = 0;
				for (int i = 0; i < junctionPositions.Count; i++)
				{
					double dist = ptvTotal.MaxDistToStructureExit(junctionPositions[i]);
					if (dist > maxDistance)
					{
						maxDistance = dist;
					}
				}

				// break if condition is satisfied or if maximum acceptable isos is reached whether the condition is met or not
				if (delta + minOverlap < 2 * fieldSizeFFS * (Machine.SID - (maxDistance)) / Machine.SID || nriso == 8)
				{
					break;
				}
				else
				{
					isos.Clear();
					junctionPositions.Clear();
				}
			}

			return isos;
		}


		// should really check the Z_PTV:s as these are probably created by using a VOI and boolean, much harder to do though.
		public void TMIcheckIfCropped(PlanSetup plan)
		{
			StructureSet ss = plan.StructureSet;
			Image image = ss.Image;

			Structure ptvTotal = ss.Structures.Where(s => s.Id == "PTV_Total").SingleOrDefault();  // TODO: better to take dicom type
			Structure body = ss.Structures.Where(s => s.Id == "BODY").SingleOrDefault();  // TODO: better to take dicom type, can only have one structure with dicom type BODY (?)
																						  //Bolus is unfortunately not a structure, can not check it the same way. However, if PTV is cropped, then it's likely the bolus is as well.

			VVector imageOrigin = image.Origin;

			double bodyLeft = body.MeshGeometry.Positions.Max(p => p.X);
			double bodyRight = body.MeshGeometry.Positions.Min(p => p.X);
			double bodyDors = body.MeshGeometry.Positions.Max(p => p.Y);
			double bodyVent = body.MeshGeometry.Positions.Min(p => p.Y);
			double bodyCran = body.MeshGeometry.Positions.Max(p => p.Z);
			double bodyCaud = body.MeshGeometry.Positions.Min(p => p.Z);

			double ptvLeft = ptvTotal.MeshGeometry.Positions.Max(p => p.X);
			double ptvRight = ptvTotal.MeshGeometry.Positions.Min(p => p.X);
			double ptvDors = ptvTotal.MeshGeometry.Positions.Max(p => p.Y);
			double ptvVent = ptvTotal.MeshGeometry.Positions.Min(p => p.Y);
			double ptvCran = ptvTotal.MeshGeometry.Positions.Max(p => p.Z);
			double ptvCaud = ptvTotal.MeshGeometry.Positions.Min(p => p.Z);

			string message = string.Empty;


			// dose matrix coordinate system is dependent on plan treatment orientation
			// can only be checked if calculated, found no way to check the calculation volume before calculation

			if (plan.Dose != null)
			{
				VVector doseOrigin = plan.Dose.Origin;
				double doseCran;
				double doseCaud;
				double doseLeft;
				double doseRight;

				if (plan.TreatmentOrientation == PatientOrientation.HeadFirstSupine)
				{
					doseCran = doseOrigin.z + plan.Dose.ZSize * plan.Dose.ZRes;
					doseCaud = doseOrigin.z;
					doseLeft = doseOrigin.x + plan.Dose.XSize * plan.Dose.XRes;
					doseRight = doseOrigin.x;
				}
				else
				{
					doseCran = doseOrigin.z;
					doseCaud = doseOrigin.z - plan.Dose.ZSize * plan.Dose.ZRes;
					doseLeft = doseOrigin.x;
					doseRight = doseOrigin.x - plan.Dose.XSize * plan.Dose.XRes;
				}

				double doseDors = doseOrigin.y + plan.Dose.YSize * plan.Dose.YRes;
				double doseVent = doseOrigin.y;

				message += "\nDose matrix margin from PTV is:\n";

				message += "Cran:\t" + (doseCran - ptvCran).ToString("0") + " mm.\n";
				message += "Caud:\t" + (-(doseCaud - ptvCaud)).ToString("0") + " mm.\n";
				message += "Left:\t" + (doseLeft - ptvLeft).ToString("0") + " mm.\n";
				message += "Right:\t" + (-(doseRight - ptvRight)).ToString("0") + " mm.\n";
				message += "Dors:\t" + (doseDors - ptvDors).ToString("0") + " mm.\n";
				message += "Vent:\t" + (-(doseVent - ptvVent)).ToString("0") + " mm.\n";
			}

			double imz = image.XSize * image.XRes;

			message += " Image size X: " + imz.ToString("0.0") + "mm\n";

			message += " Distance from BODY to PTV_Total (maximum position in each direction): \n";

			if (plan.TreatmentOrientation == PatientOrientation.HeadFirstSupine)
			{
				message += "Cran:\t" + (ptvCran - bodyCran).ToString("0") + " mm.\n";
			}
			else
			{
				message += "Caud:\t" + (-(ptvCaud - bodyCaud)).ToString("0") + " mm.\n";
			}

			message += "Left:\t" + (ptvLeft - bodyLeft).ToString("0") + " mm.\n";
			message += "Right:\t" + (-(ptvRight - bodyRight)).ToString("0") + " mm.\n";
			message += "Dors:\t" + (ptvDors - bodyDors).ToString("0") + " mm.\n";
			message += "Vent:\t" + (-(ptvVent - bodyVent)).ToString("0") + " mm.\n";

			MessageBox.Show(message);

		}

		/// <summary>
		/// Determine the isocenter position necessary to cover the cranial PTV in all angles for 90 deg coll and x2 13.5 cm
		/// </summary>
		/// <param name="plan"></param>
		/// <returns></returns>
		private string TMICranialCoverage(PlanSetup plan)
		{
			// TODO: need to check if all iso are positioned in same lat and vrt value first
			StructureSet ss = plan.StructureSet;
			Image image = ss.Image;
			Structure structure = ss.Structures.Where(s => s.Id == "PTV_Total").SingleOrDefault();  // TODO: better to take dicom type, might want to use PTV_Total instead
			string result = string.Empty;
			Beam cranBeam = plan.Beams.Where(b => b.Id == "1").FirstOrDefault();
			VVector isoPos = cranBeam.IsocenterPosition;
			double assumedFieldsize = 135;  // TODO: Check the actual field size in cranial direction...    **********************************double x1 = cranBeam.ControlPoints.First. ... nope

			string debug = string.Empty;

			// check the max distance from junction/iso to the four corners of the bounds of the structure, same values for all junctions
			double searchRadius = structure.MaxDistToBoundsExit(isoPos);

			// Rotate around the structure with origo in iso and check intersection with structure. Report optimal position
			// of iso to cover PTV. Enough to scan from iso to outer position determined by scan radius, divergence and field size

			int anglesIntersecting = 0;

			double zOffsetIso = assumedFieldsize;
			double zOffsetEnd = assumedFieldsize * (Machine.SID - searchRadius) / Machine.SID;

			anglesIntersecting = NrAnglesIntersectStructure(structure, isoPos, searchRadius, zOffsetIso, zOffsetEnd);



			switch (anglesIntersecting)
			{
				case 0:
					result += "Cranial part of PTV covered ok.\n";
					break;
				case 360:
					result += "** Cranial part of PTV not covered!\n";
					break;
				default:
					result += "Cranial part of PTV not covered in all angles." + "(" + anglesIntersecting + "/360)\n";
					break;
			}

			VVector optimalPos = image.DicomToUser(ZposCoverAllAngles(structure, isoPos, assumedFieldsize, true), plan);
			result += "To precisely cover the cranial part of PTV_Total in all angles, given the planned isocenter position in lat and vrt, the z-position is: ";
			result += (optimalPos.z / 10).ToString("0.0") + " cm\n";
			return result;
		}


		private string TMICheckForCouchValues(PlanSetup plan)
		{
			string cResults = string.Empty;
			string beamIdWithCouchValues = string.Empty;
			foreach (var beam in plan.Beams)
			{
				if (!beam.ControlPoints[0].TableTopLateralPosition.ToString().Contains("NaN") || !beam.ControlPoints[0].TableTopLongitudinalPosition.ToString().Contains("NaN") || !beam.ControlPoints[0].TableTopVerticalPosition.ToString().Contains("NaN"))
				{
					beamIdWithCouchValues += beam.Id + ", ";
				}
			}

			if (beamIdWithCouchValues.Length > 0)
			{
				cResults = "* Please remove couch coordinates from the following fields: " + beamIdWithCouchValues.Substring(0, beamIdWithCouchValues.Length - 2) + ".\n\n";
			}
			return cResults;
		}



		/// <summary>
		/// returns z-position of isocenter in dicom coordinates where given field size [mm] covers the structure in all angles
		/// either cranially or caudally
		/// </summary>
		/// <param name="structure"></param>
		/// <param name="iso"></param>
		/// <param name="fieldSize"></param>
		/// <returns></returns>
		public VVector ZposCoverAllAngles(Structure structure, VVector iso, double fieldSize, bool cranial)
		{
			double searchRadius = structure.MaxDistToBoundsExit(iso);
			double zOffsetIso = fieldSize;
			double zOffsetEnd = fieldSize * (Machine.SID - searchRadius) / Machine.SID;
			int startScanInStructure = 5; // offset value for start scanning; mm before structure end in long
			int maxScanDist = 200; // high value needed if offset in transverce plane large
			VVector optimalPosDicom = iso;
			double anglesIntersecting;
			if (cranial)
			{
				optimalPosDicom.z = structure.MeshGeometry.Positions.Max(p => p.Z) - startScanInStructure - fieldSize;
				for (int i = 0; i < maxScanDist; i++)
				{
					anglesIntersecting = NrAnglesIntersectStructure(structure, optimalPosDicom, searchRadius, zOffsetIso + i, zOffsetEnd + i);
					if (anglesIntersecting == 0)
					{
						optimalPosDicom.z += i;
						break;
					}
				}
			}
			else // caudal
			{
				optimalPosDicom.z = structure.MeshGeometry.Positions.Min(p => p.Z) + startScanInStructure + fieldSize;
				for (int i = 0; i < maxScanDist; i++)
				{
					anglesIntersecting = NrAnglesIntersectStructure(structure, optimalPosDicom, searchRadius, -zOffsetIso - i, -zOffsetEnd - i);
					if (anglesIntersecting == 0)
					{
						optimalPosDicom.z -= i;
						break;
					}
				}
			}
			return optimalPosDicom;
		}


		/// <summary>
		/// Calculates nr of angles that intersects a structure when profiling from iso (+offset z) to endpoint defined by
		/// searchradius and offset z from iso. Parameters in mm and Dicom coordinates.
		/// </summary>
		/// <param name="structure"></param>
		/// <param name="beam"></param>
		/// <param name="searchRadius"></param>
		/// <param name="zOffsetStart"></param>
		/// <param name="zOffsetEnd"></param>
		/// <param name="angleIncrement"></param>
		/// <returns></returns>
		public int NrAnglesIntersectStructure(Structure structure, VVector iso, double searchRadius, double zOffsetStart, double zOffsetEnd, double angleIncrement = 1)
		{
			VVector startPoint = iso;
			VVector endPoint = iso;
			startPoint.z += zOffsetStart;
			endPoint.z += zOffsetEnd;

			int nrAnglesIntersecting = 0;
			for (double angle = 0; angle < 360; angle += angleIncrement)
			{
				endPoint.x = iso.x + searchRadius * Math.Cos(angle * Math.PI / 180);
				endPoint.y = iso.y + searchRadius * Math.Sin(angle * Math.PI / 180);

				var structureOutside = structure.GetSegmentProfile(startPoint, endPoint, new BitArray(100)).Where(s => s.Value == false);
				if (structureOutside.Count() != 100)
				{
					nrAnglesIntersecting++;
				}
			}
			return nrAnglesIntersecting;
		}


		/// <summary>
		/// Returns a list of VVectors representing all the isocenters found in plan ordered by iso position z, then x, then y
		/// </summary>
		/// <param name="plan"></param>
		/// <returns></returns>
		public List<VVector> TMIIsocentersInPlan(PlanSetup plan)
		{
			List<VVector> isos = new List<VVector>();
			foreach (var beam in plan.Beams.OrderByDescending(b => b.IsocenterPosition.z).ThenBy(b => b.IsocenterPosition.x).ThenBy(b => b.IsocenterPosition.y))
			{
				if (isos.Any(i => i.x == beam.IsocenterPosition.x && i.y == beam.IsocenterPosition.y && i.z == beam.IsocenterPosition.z) == false)
				{
					isos.Add(beam.IsocenterPosition);
				}
			}
			return isos;
		}


		private string CheckBodyCenter(PlanSetup plan)
		{
			string results = string.Empty;

			StructureSet ss = plan.StructureSet;
			Image image = ss.Image;
			Structure body = ss.Structures.Where(s => s.Id == "BODY").SingleOrDefault();  // TODO: better to take dicom type
			if (body.MeshGeometry == null)
			{
				results = "Cannot find geometric center";
			}
			else
			{
				VVector gCenter = body.GeometricCenter();
				VVector gCEclipse = image.DicomToUser(gCenter, plan);
				VVector wCenter = body.CenterPoint;
				VVector wCEclipse = image.DicomToUser(wCenter, plan);
				results = "Body weighted center: \tvrt: " + (wCEclipse.y / 10).ToString("0.0") + "\tlat: " + (wCEclipse.x / 10).ToString("0.0") + "\n" +
						"Body geometric center: \tvrt: " + (gCEclipse.y / 10).ToString("0.0") + "\tlat: " + (gCEclipse.x / 10).ToString("0.0") + "\n";
			}

			return results;
		}


		/// <summary>
		/// Information for sum plan exclusively for TBI/TMI plans with delta shifts between plan isocenters
		/// </summary>
		/// <param name="psum"></param>
		private void CheckPlanSum(PlanSum psum)
		{// Only for TMI...
			string sumPlans = string.Empty;
			PlanSetup firstPlan = psum.PlanSetups.OrderBy(p => p.Id).FirstOrDefault();
			PlanCat planCat = PlanCat.Unknown;
			CheckPlanCategory(firstPlan, ref planCat);
			if (planCat == PlanCat.TBI || planCat == PlanCat.TMI)
			{
				//List<PlanSetup> sumPlanTreatOrderFFS = psum.PlanSetups.Where(p => p.TreatmentOrientation == PatientOrientation.FeetFirstSupine).OrderBy(p => p.Beams.Where(b => b.IsSetupField).Count()).ThenBy(p => p.Id).ToList();
				List<PlanSetup> sumPlanTreatOrderHFS = psum.PlanSetups.Where(p => p.TreatmentOrientation == PatientOrientation.HeadFirstSupine).OrderBy(p => p.Id).ToList();
				List<PlanSetup> sumPlanTreatOrderFFS = psum.PlanSetups.Where(p => p.TreatmentOrientation == PatientOrientation.FeetFirstSupine).OrderBy(p => p.Id).ToList();

				// Need to check that it actually is the treatment plan, how ? name of the structure set CX DP
				if (sumPlanTreatOrderHFS.Count() >= 1)
				{
					sumPlans += "HFS" + "\n\n" + sumPlanTreatOrderHFS[0].Id + "\t\n" + DeltaShiftFromOrigin(sumPlanTreatOrderHFS[0]);
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
			}
			else
			{
				sumPlans = "Can only show sumplans for category TMI and TBI VMAT";
			}


			MessageBox.Show(sumPlans);
		}
		#endregion

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
				return PlanCat.TBI; // Oldschool TBI
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
			if (Math.Round(maxRadi, 0) == Math.Round(minRadi, 0) && standardEBlockDiameters.Contains((int)(2 * Math.Round(maxRadi, 0) / 10)))
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



		#region delta shift

		private string DeltaShiftFromPlanToPlan(PlanSetup fromPlan, PlanSetup toPlan)
		{
			VVector deltaShift = DeltaShiftIncm(fromPlan, fromPlan.Beams.First().IsocenterPosition, toPlan.Beams.First().IsocenterPosition);

			string delta = " Delta(Vrt,Lng,Lat)[cm]: \t" + deltaShift.y.ToString("0.00");
			delta += "\t" + deltaShift.z.ToString("0.00");
			delta += "\t" + deltaShift.x.ToString("0.00") + "\n";
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
		private VVector DeltaShiftOrigin(PlanSetup plan)
		{
			return DeltaShiftIncm(plan, plan.StructureSet.Image.UserOrigin, plan.Beams.First().IsocenterPosition);
		}
		private VVector DeltaShiftPlanToPlan(PlanSetup fromPlan, PlanSetup toPlan)
		{
			return DeltaShiftIncm(fromPlan, fromPlan.Beams.First().IsocenterPosition, toPlan.Beams.First().IsocenterPosition);
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

		#region Bolus // thickness in ID still used in 16.1


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
				if (Double.TryParse(g.Value.Replace(",", "."), out bolusThick))
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


			CheckStructNameAndType(plan);



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
			List<Structure> ptvList = ss.Structures.Where(s => s.Id.Substring(0, 3) == "PTV").Where(s => s.DicomType == "PTV").ToList();
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



		private void CheckStructNameAndType(PlanSetup plan)
		{
			// From dose value and plan dose recommend isodose lines?

			// SSM 2016:18 national nomenclature for structures in radiation therapy (ignores the fact that "/" also is ok for "free text" at the end)
			Regex targetNaming = new Regex(@"^([GCIP]TV)(?<type>[TNM]?)(?<number>[1-9]?)(_(?<position>[LR]))?(_(?<dose>\d{1,2}(?:\.[0-9]{1,2})?))?(?:[\(](?<description>.+)[\)])?", RegexOptions.IgnoreCase);


			/*
			 * Order: by Id:  PTV, E_PTV, ITV, CTV, GTV, BODY, BolusIncluded in calc,	 avoidance, organs and the rest, help structures, Bolus not included,	helpBolus, empty structures
			 * 
			 * Dicom type	description		volume		nrofparts	assigned hu		comment in eclipse	
			 */




			StructureSet ss = plan.StructureSet;
			// structures included resp excluded in calculation (outside body and not bolus, empty)
			// sort order type; ptv, ctv, itv, gtv 
			List<Structure> structures = new List<Structure>();
			List<Structure> EmptyStr = new List<Structure>();

			structures = ss.Structures.Where(s => s.HasSegment).OrderBy(s => s.Id).ToList();
			EmptyStr = ss.Structures.Where(s => s.HasSegment == false).ToList();

			//int nrOfPTV = ss.Structures.Where(s => s.Id.Length > 4).Where(s => s.Id.Substring(0, 3) == "PTV").Where(s => s.DicomType == "PTV").Count();
			//structures.OrderByDescending(s => s.Id.Substring(0, 3) == "PTV").ThenByDescending(s => s.Id.Substring(0, 3) == "CTV").ThenByDescending(s => s.Id.Substring(0, 3) == "GTV").ToList();

			// check dicom types and structure codes for the structures named as target volumes

			List<Structure> ptvStruct = structures.Where(s => s.Id.Length >= 3 && s.Id.Substring(0, 3).ToUpper() == "PTV").OrderBy(s => s.Id).ToList();
			List<Structure> ctvStruct = structures.Where(s => s.Id.Length >= 3 && s.Id.Substring(0, 3).ToUpper() == "CTV").OrderBy(s => s.Id).ToList();
			List<Structure> itvStruct = structures.Where(s => s.Id.Length >= 3 && s.Id.Substring(0, 3).ToUpper() == "ITV").OrderBy(s => s.Id).ToList();
			List<Structure> gtvStruct = structures.Where(s => s.Id.Length >= 3 && s.Id.Substring(0, 3).ToUpper() == "GTV").OrderBy(s => s.Id).ToList();

			List<Structure> targetStructures = new List<Structure>();
			targetStructures.AddRange(ptvStruct);
			targetStructures.AddRange(ctvStruct);
			targetStructures.AddRange(itvStruct);
			targetStructures.AddRange(gtvStruct);



			string message = string.Empty;
			foreach (var target in targetStructures)
			{
				string searchString = target.Id.Substring(0, 3).ToUpper();
				if (target.DicomType.Contains(searchString) == false)
				{
					message += target.Id + ": ändra till volymtyp " + searchString + " (är nu satt som " + target.DicomType + ")\n";
				}
				if (target.StructureCode.Code.Contains(searchString) == false)
				{
					message += target.Id + ": ändra till strukturkod " + searchString + " (är nu satt som " + target.StructureCode.Code.FirstOrDefault() + ")\n";
				}
			}

			//*************  check of naming against SSM 2016:18 national standards

			// cant get it to work solely on regex for som reason, need to compare string lengths as well
			foreach (var target in targetStructures)
			{
				if (targetNaming.IsMatch(target.Id) && targetNaming.Match(target.Id).Length == target.Id.Length)
				{
					continue;
				}
				else
				{
					message += target.Id + ": namgivning verkar inte följa nationell standard. \n";
				}
			}

			Structure body = ss.Structures.Where(s => s.DicomType == "EXTERNAL").SingleOrDefault();  // ***************** crashes if no BODY type found!
			VVector bodyCenter = body.GeometricCenter();
            // TODO: generalize this for all structures. check position of target(L/ R) compared to body center

            foreach (var target in targetStructures)
            {
                if (targetNaming.Match(target.Id).Length == target.Id.Length && !string.IsNullOrEmpty(targetNaming.Match(target.Id).Groups["position"].Value))
                {

                    switch (targetNaming.Match(target.Id).Groups["position"].Value)
                    {

                        case "L":

                            if (target.GeometricCenter().x < bodyCenter.x)
                            {
                                message += target.Id + "Kontrollera namngivning, geometriskt center verkar vara till höger om BODY-center\n";
                            }

                            break;
                        case "R":

                            if (target.GeometricCenter().x > bodyCenter.x)
                            {
                                message += target.Id + "Kontrollera namngivning, geometriskt center verkar vara till vänster om BODY-center\n";
                            }

                            break;
                        default:
                            break;
                    }
                }
            }


            // ******************************   testing testing *********************

            foreach (var ptv in ptvStruct)
			{
				double doseValue = 0.0;
				if (targetNaming.IsMatch(ptv.Id) && string.IsNullOrEmpty(targetNaming.Match(ptv.Id).Groups["dose"].Value) == false && targetNaming.Match(ptv.Id).Length == ptv.Id.Length)
				{
					string doseValueString = targetNaming.Match(ptv.Id).Groups["dose"].Value;
					if (double.TryParse(doseValueString, out doseValue))
					{
						message += ptv.Id + " testing: dose value = " + doseValue.ToString("0.00") + " Gy.\n\n";
					}
				}
				else
				{
					message += ptv.Id + " ingen dosnivå i PTV-ID.\n\n";
				}

				double caud = ptv.MeshGeometry.Bounds.Z;
				message += "Testing: " + ptv.Id + ": starts from slice " + ss.Image.ClosestSlice(caud) + "\n";
			}












			if (string.IsNullOrEmpty(message) == false)
			{
				MessageBox.Show(message);
			}



			//DICOM types
			//Possible values are "AVOIDANCE", "CAVITY", "CONTRAST_AGENT", "CTV", "EXTERNAL", "GTV", "IRRAD_VOLUME", 
			//"ORGAN", "PTV", "TREATED_VOLUME", "SUPPORT", "FIXATION", "CONTROL", and "DOSE_REGION". 


			string listOfStr = string.Empty;
			foreach (var str in structures)
			{
				listOfStr += str.Id + "\t";//						 + str.DicomType + "\t" + str.StructureCode + "\n";
			}

			//MessageBox.Show(listOfStr);
		}











		//TODO: create method for getting coord 8 corners of max pos CVT  and check if the coord is included within respective 8 corners of PTV
		private string SBRTCtvPtvPlanNumbering(PlanSetup plan, PlanCat planCat)
		{
			string cResults = string.Empty;
			StructureSet ss = plan.StructureSet;
			int number;

			// check that CTV, PTV and Plan numbering is consistent, only necessary if more than one PTV
			// PTV from plans target volume, ctv from PTV number and check that mass center is within respective PTV boundaries
			// Plan numbering not necessarily consistent with PTV if multiple structure sets are used in the same course 
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

		/// <summary>
		/// Special check for SBRT structure set to, if necessary, remind the user to include fixation devices in the BODY structure,
		/// e.g. vacuum bag, SBF etc. Also checks the slice thickness against guidelines.
		/// </summary>
		/// <param name="plan"></param>
		/// <returns></returns>
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
				Structure body = ss.Structures.Where(s => s.Id == "BODY").SingleOrDefault();  // TODO: better to take dicom type
				double bodyVrtMin = body.MeshGeometry.Bounds.Y + body.MeshGeometry.Bounds.SizeY;
				int minDifferenceBodySkin = 3; // if difference larger than 3 mm assume that body includes fixation device

				VVector endPoint = plan.Beams.FirstOrDefault().IsocenterPosition;
				VVector startPoint = endPoint;
				startPoint.y = bodyVrtMin;

				var bodyEntrance = body.GetSegmentProfile(startPoint, endPoint, new BitArray(100)).Where(x => x.Value == true).First();
				var skinEntrance = skin.GetSegmentProfile(startPoint, endPoint, new BitArray(100)).Where(x => x.Value == true).First();

				if (bodyEntrance.Position.y <= skinEntrance.Position.y + minDifferenceBodySkin)
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
		//TODO: first find structure in structure set whitch corresponds to target structure, then do check on that
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

					target = sSet.Structures.Where(s => s.Id == plan.TargetVolumeID).SingleOrDefault();
					cResults = "* Plan target volume should be of type PTV.\nPlan target volume: " + target.Id + " of Dicom-type: " + target.DicomType + ".\n";
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

			if (char.ToUpperInvariant(planIdFirstChar) == 'P' && int.TryParse(planIdSecondChar.ToString(), out planNumber))
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
		// Om enbart CBCT-fält, beroende på plancat, kan ge hint om ext eller flytt från cbctmatchstruktur, plan target volume, iso 
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
				if (beam.Id.Length >= 2   &&     beam.Id.ToUpper().Substring(0, 2).Equals(plan.Id.ToUpper().Substring(0, 2))) // naming convention, should start with first two char in plan ID
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
				else if (plan.ApprovalStatus == PlanSetupApprovalStatus.PlanningApproved || plan.ApprovalStatus == PlanSetupApprovalStatus.TreatmentApproved)
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

		private string CheckFieldRules(PlanSetup plan, PlanCat planCat)
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
				remarks += TMICheckForCouchValues(plan);
			}




			foreach (var beam in plan.Beams.Where(b => !b.IsSetupField).OrderBy(b => b.Id))
			{
				//BeamExtraInfo beamExtras = new BeamExtraInfo(beam);

				cResults = cResults + beam.Id.PadRight(2) + "\t" + beam.EnergyModeDisplayName.PadRight(6) + "\t" + beam.Technique.Id + "\t" + beam.DoseRate.ToString().PadRight(4) + "\t" + beam.MLCPlanType + "\t";
				//cResults += Math.Round(BeamExtensions.EstimatedBeamOnTime(beam)).ToString("0") + " s\t" + (beam.ControlPoints.Count()) +"\n";
				cResults += beam.EstimatedBeamOnTime() + " s\n";
				beamOnTimeInSec += BeamExtensions.EstimatedBeamOnTime(beam);



				if (planCat == PlanCat.SBRT && beam.Technique.Id.Contains("SRS") == false)
				{
					if (countSRSRemarks < 1) { remarks = remarks + "** Change technique to SRS-" + beam.Technique.Id + "! \n"; };
					countSRSRemarks++;
				}
				if (beam.EnergyModeDisplayName.Contains("FFF"))
				{
					remarks += CheckDoseRateFFF(beam, ref countDoseRateFFFRemarks);
				}
				if (beam.MLCPlanType == MLCPlanType.ArcDynamic && planCat == PlanCat.SBRT)
				{
					remarks += CheckArcDynFFF(beam, ref countArcDynFFFRemarks);
					remarks += SBRTCheckArcDynCollAngle(beam, ref countArcDynCollAngleRemarks);
				}
				if (beam.MLCPlanType == MLCPlanType.Static)
				{
					if (planCat == PlanCat.SBRT)
					{
						remarks += CheckMLCStaticFFF(plan, beam, ref countMLCStaticFFFRemarks);
					}
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
			return cResults + "\n" + "Estimated total beam-on-time: " + (beamOnTimeInSec / 60).ToString("0.0") + " min\n\n" + "Estimated delivery time: > " + Math.Round(GetEstimatedTreatmentTime(plan) / 60, 1).ToString("0.0") + " min\n" + remarks;
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
				treatTime += BeamExtensions.EstimatedBeamOnTime(beam);
			}

			// control sequens and energy change happens parallel with, and independent to, mechanical movements
			for (int i = 1; i < cpFirst.Count(); i++)
			{
				mechMovementTime.Add(DeltaAngle(cpLast[i - 1].GantryAngle, cpFirst[i].GantryAngle) / Machine.GantryMaxSpeed);
				mechMovementTime.Add(MaxDeltaJaw(cpLast[i - 1], cpFirst[i]) / Machine.JawXMaxSpeed);
				if (beamsInOrder[i - 1].EnergyModeDisplayName.Equals(beamsInOrder[i]))
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

				if (plan.AutomationPrerequisites() == false)
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


		// ********* 	Kontroll av dosrat vid FFF	*********

		public string CheckDoseRateFFF(Beam beam, ref int countDoseRateFFFRemarks)
		{
			string cResults = "";

			if (beam.MaxDoseRateUsed() == false)
			{
				if (countDoseRateFFFRemarks < 1)
				{
					cResults = "* Consider changing dose rate to maximum.\n";
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
			if (countMLCStaticFFFRemarks < 1 && !beam.EnergyModeDisplayName.Contains("FFF") && plan.AnyWedges() == false)
			{
				cResults = "* Consider changing energy to FFF \n";
				countMLCStaticFFFRemarks++;
			}
			return cResults;
		}



		// ********* 	Kontroll av kollimatorvinkel vid Dynamic Arc	*********

		public string SBRTCheckArcDynCollAngle(Beam beam, ref int countArcDynCollAngleRemarks)
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
