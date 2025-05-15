//
//
// using namespace adsk::core;
// using namespace adsk::fusion;
//
// // ------------------------------------------------
// // external helper supplied by the user
// extern std::string saveAirfoil(const std::string& nacaCode);
//
// // ----- USER INPUT DATA -------------------------------------------------
// static const std::vector<std::string>  nacaCodes = { "2412", "2412", "0012",
//                                                      "0010", "0008", "0006" };
// static const std::vector<double> spanPos_mm = { 0, 200, 400, 600, 800, 1000 }; // mm
// static const std::vector<double> chord_mm   = { 150, 145, 135, 120, 100,  60 }; // mm
// static const std::vector<double> twist_deg  = {   0,   -2,  -4,  -6,  -8, -10 }; // °
// static const bool keepDesignOpen = true; // set false to close design on finish
// // -----------------------------------------------------------------------
//
// static Ptr<Application>  app;
// static Ptr<UserInterface> ui;
//
// struct Point2D { double x, y; };
//
// // tiny 2×2 rotation
// static void rotateYz(double& y, double& z, double angleRad)
// {
//     double c = std::cos(angleRad), s = std::sin(angleRad);
//     double ty =  c*y + s*z;
//     double tz = -s*y + c*z;
//     y = ty; z = tz;
// }
//
// // ----------------------------------------------------------
// // Parse a simple SECTION-less *.dat/*.txt file.
// // Assumes x y pairs, first = last for closed loop.
// static std::vector<Point2D> readAirfoil(const std::string& path)
// {
//     std::ifstream f(path);
//     std::vector<Point2D> pts;
//     if(!f) return pts;
//     std::string line;
//     while(std::getline(f, line))
//     {
//         if(line.empty() || line[0]=='#') continue;
//         std::istringstream iss(line);
//         Point2D p;
//         if(iss >> p.x >> p.y) pts.push_back(p);
//     }
//     return pts;
// }
//
// // ----------------------------------------------------------
// // Create one sketch containing a closed fitted spline
// // and return its profile.
// static Ptr<Profile> createSectionSketch(Ptr<Component> root,
//                                         double offsetXmm,
//                                         const std::vector<Point2D>& section,
//                                         double chord,
//                                         double twistDeg)
// {
//     if(section.size() < 3) return nullptr;
//
//     // ------------------------------------------------ plane at span position
//     Ptr<ConstructionPlanes> planes = root->constructionPlanes();
//     Ptr<ConstructionPlane>  baseXY = root->xYConstructionPlane();
//     Ptr<ConstructionPlaneInput> pIn = planes->createInput();
//     pIn->setByOffset(baseXY,
//         ValueInput::createByReal(offsetXmm/* cm in API units = mm here /10 */ / 10.0));
//     Ptr<ConstructionPlane> plane = planes->add(pIn);
//
//     // ------------------------------------------------ new sketch
//     Ptr<Sketches> sketches = root->sketches();
//     Ptr<Sketch> sk = sketches->add(plane);
//
//     // ------------------------------------------------ build ObjectCollection
//     Ptr<ObjectCollection> oc = ObjectCollection::create();
//     double twistRad = twistDeg * M_PI/180.0;
//
//     for(const auto& p : section)
//     {
//         // scale and twist
//         double y =  p.x * chord;        // chordwise axis -> Y (mm)
//         double z =  p.y * chord;        // thickness axis -> Z (mm)
//         rotateYz(y, z, twistRad);
//
//         Ptr<Point3D> pt = Point3D::create(offsetXmm/10.0, y/10.0, z/10.0); // cm
//         oc->add(pt);
//     }
//     // fitted spline (closed = first==last)
//     Ptr<SketchFittedSpline> spline =
//         sk->sketchCurves()->sketchFittedSplines()->add(oc);
//
//     if(!spline) return nullptr;
//     return sk->profiles()->item(0);
// }
//
// // ----------------------------------------------------------
// // Main routine executed by Fusion when the add-in is run
// extern "C" XI_EXPORT bool run(const char* /*context*/)
// {
//     app = Application::get();
//     if(!app) return false;
//     ui  = app->userInterface();
//
//     Ptr<Design> design = app->activeProduct();
//     if(!design)
//     {
//         ui->messageBox("No active design");
//         return false;
//     }
//     Ptr<Component> root = design->rootComponent();
//
//     // ------------------------------------------------ generate all sections
//     std::vector< Ptr<Profile> > loftProfiles;
//     for(size_t i=0; i<nacaCodes.size(); ++i)
//     {
//         std::string file = saveAirfoil(nacaCodes[i]); // <== user-supplied
//         auto sectionPts  = readAirfoil(file);
//         auto profile = createSectionSketch(root,
//                                            spanPos_mm[i],
//                                            sectionPts,
//                                            chord_mm[i],
//                                            twist_deg[i]);
//         if(profile)
//             loftProfiles.push_back(profile);
//         else
//             ui->messageBox("Section "+std::to_string(i)+" failed.");
//     }
//
//     // ------------------------------------------------ create loft
//     if(loftProfiles.size() >= 2)
//     {
//         Ptr<ObjectCollection> loftSecs = ObjectCollection::create();
//         for(auto& p : loftProfiles) loftSecs->add(p);
//
//         Ptr<LoftFeatures> lofts = root->features()->loftFeatures();
//         Ptr<LoftFeatureInput> lIn = lofts->createInput(
//             loftSecs, FeatureOperations::NewBodyFeatureOperation);
//         lIn->isSolid(true);
//         lofts->add(lIn);
//     }
//     else
//         ui->messageBox("Need at least two valid sections to loft.");
//
//     if(!keepDesignOpen) app->closeDocument(design->parentDocument());
//
//     return true;   // success
// }
//
// // ----------------------------------------------------------
// extern "C" XI_EXPORT bool stop(const char* /*context*/)
// {
//     return true;
// }