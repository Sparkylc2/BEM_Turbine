//
// Created by Lukas Campbell on 14/05/2025.
//

#pragma once


using namespace SldWorks;
using namespace std::string_literals;

// ---------------------------------------------------------------------------
// Utility: read 2-column coordinate file => vector<double x, double y>
std::vector<std::pair<double, double>>
loadPoints(const std::filesystem::path& file)
{
    std::ifstream in(file);
    if (!in) throw std::runtime_error("Cannot open " + file.string());

    std::vector<std::pair<double, double>> pts;
    double x{}, y{};
    while (in >> x >> y) pts.emplace_back(x, y);
    return pts;
}

// ---------------------------------------------------------------------------
// Sketch the closed polyline on the current sketch plane
void drawAirfoilSketch(ISketchManager* smgr,
                       const std::vector<std::pair<double,double>>& pts,
                       double chord, double twist_deg)
{
    const double deg2rad = M_PI / 180.0;
    double cosT = std::cos(twist_deg * deg2rad);
    double sinT = std::sin(twist_deg * deg2rad);

    // Begin & end sketch editing
    CComPtr<ISketch> sketch;
    smgr->InsertSketch(true);                  // true = add new
    smgr->get_ActiveSketch(&sketch);

    long n = static_cast<long>(pts.size());
    for (long i = 0; i < n; ++i) {
        auto [x0,y0] = pts[i];
        auto [x1,y1] = pts[(i+1)%n];           // wrap to first point

        // scale & twist
        double x0c =  chord*( x0*cosT - y0*sinT );
        double y0c =  chord*( x0*sinT + y0*cosT );
        double x1c =  chord*( x1*cosT - y1*sinT );
        double y1c =  chord*( x1*sinT + y1*cosT );

        VARIANT vLine;
        smgr->CreateLine(x0c, y0c, 0.0, x1c, y1c, 0.0, &vLine);
    }

    smgr->InsertSketch(true);                  // exit sketch
}

// ---------------------------------------------------------------------------
// MAIN
// ---------------------------------------------------------------------------
int wmain()
{
    // ------------------------------------------------------------------ user data
    const std::vector<std::string> nacaCodes = {
        "2412","4415","0012","6412","2424","9412","3410"
    };
    std::vector<double> spanPos  = {0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0}; // mm
    std::vector<double> chords   = {1.00, 0.95, 0.90, 0.80, 0.70, 0.60, 0.50};         // scale
    std::vector<double> twists   = { 0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0};          // deg
    // -----------------------------------------------------------------------------

    if (nacaCodes.size()!=spanPos.size() ||
        spanPos.size()!=chords.size()     ||
        chords.size()!=twists.size())
        throw std::runtime_error("Vector length mismatch!");

    // Initialise COM & launch / connect to SolidWorks
    CoInitialize(nullptr);
    {
        CComPtr<ISldWorks> swApp;
        HRESULT hr = swApp.CoCreateInstance(__uuidof(SldWorks::SldWorks), nullptr,
                                            CLSCTX_LOCAL_SERVER);
        if (FAILED(hr) || !swApp) throw std::runtime_error("Cannot start SolidWorks");

        swApp->Visible = VARIANT_TRUE;

        // Create a new part document
        CComBSTR tmpl(L"");   // empty => default template
        CComPtr<IModelDoc2> model;
        swApp->NewDocument(tmpl, 0, 0.0, 0.0, &model);

        CComPtr<IFeatureManager> featMgr;
        model->get_FeatureManager(&featMgr);

        // For collecting section sketches
        std::vector<CComPtr<IRefPlane>> planes;
        std::vector<CComPtr<IEntity>>    profiles;

        // Loop over every airfoil section
        for (size_t i = 0; i < nacaCodes.size(); ++i) {
            // 1) Produce & load coordinates ----------------------------------
            std::filesystem::path file = saveAirfoil(nacaCodes[i]);  // ← your function
            auto pts = loadPoints(file);

            // 2) Create offset plane ----------------------------------------
            CComPtr<IRefPlane> plane;
            {
                model->ClearSelection2(VARIANT_TRUE); // deselect
                // Select Front Plane (name depends on locale – safer: use GetFrontPlane)
                CComPtr<IModelDocExtension> ext;
                model->get_Extension(&ext);
                CComPtr<IEntity> frontPlaneEnt;
                ext->SelectByID2(L"Front Plane", L"PLANE", 0,0,0, VARIANT_FALSE,
                                 0, nullptr, 0, &frontPlaneEnt);
                CComVariant vDist(spanPos[i]/1000.0); // SolidWorks internal units: meters
                featMgr->InsertRefPlane((long)swRefPlaneReferenceConstraints::swRefPlaneReferenceConstraint_Distance,
                                        0, &vDist, 0, 0, 0, VARIANT_FALSE, (IDispatch**)&plane);
            }

            // 3) Draw sketch on that plane ----------------------------------
            // Select plane for sketching
            plane->Select(VARIANT_TRUE);
            CComPtr<ISketchManager> skMgr;
            model->get_SketchManager(&skMgr);
            drawAirfoilSketch(skMgr, pts, chords[i], twists[i]);

            // 4) Store sketch profile for loft ------------------------------
            CComPtr<ISelectionMgr> selMgr;
            model->get_SelectionManager(&selMgr);
            CComPtr<IEntity> skEnt;
            selMgr->GetSelectedObject6(1, -1, (IDispatch**)&skEnt); // latest sketch
            profiles.push_back(skEnt);
        }

        // --------------------------------------------------------------------
        // Create loft feature through collected profiles
        model->ClearSelection2(VARIANT_TRUE);
        for (auto& ent : profiles)
            ent->Select(VARIANT_TRUE);

        VARIANT_BOOL boolStatus{};
        long loftType = (long)swFeatureOperations_e::swFeatureOperation_Cut; // or _Protrusion
        CComPtr<IFeature> loft;
        featMgr->InsertProtrusionBlend(loftType, VARIANT_TRUE, VARIANT_TRUE,
                                       0, 0, 0, 0,
                                       VARIANT_FALSE, VARIANT_FALSE,
                                       0, 0, 0, nullptr, nullptr, &loft, &boolStatus);

        // --------------------------------------------------------------------
        // Save & leave SolidWorks open
        CComBSTR saveAs(L"lofted_wing.SLDPRT");
        VARIANT_LONG errors{}, warnings{};
        model->SaveAs3(saveAs, 0, swSaveAsOptions_e::swSaveAsOptions_Silent, &errors, &warnings);
    }
    CoUninitialize();
    return 0;
}

#endif //SOLIDWORKS_HOOK_H
