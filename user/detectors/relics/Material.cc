#include "Material.hh"

#include <G4Version.hh>
#include <G4UIcommand.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4UnitsTable.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalBorderSurface.hh>

#include <vector>

using std::vector;

MaterialRegister<Material> Material::reg("Material");

void Material::defineMaterials(const BambooParameters &parameters) {
    using namespace CLHEP;
    auto *pNistManager = G4NistManager::Instance();

    // ==== Select elements ====
    vector<G4Element *> elements;
    vector<G4Material *> materials;
    auto H = pNistManager->FindOrBuildElement(1);
    elements.push_back(H);
    auto B = pNistManager->FindOrBuildElement(5);
    elements.push_back(B);
    auto C = pNistManager->FindOrBuildElement(6);
    elements.push_back(C);
    auto N = pNistManager->FindOrBuildElement(7);
    elements.push_back(N);
    auto O = pNistManager->FindOrBuildElement(8);
    elements.push_back(O);
    auto F = pNistManager->FindOrBuildElement(9);
    elements.push_back(F);
    auto Al = pNistManager->FindOrBuildElement(13);
    elements.push_back(Al);
    auto Si = pNistManager->FindOrBuildElement(14);
    elements.push_back(Si);
    auto Ar = pNistManager->FindOrBuildElement(18);
    elements.push_back(Ar);
    auto Ca = pNistManager->FindOrBuildElement(20);
    elements.push_back(Ca);
    auto Cr = pNistManager->FindOrBuildElement(24);
    elements.push_back(Cr);
    auto Mn = pNistManager->FindOrBuildElement(25);
    elements.push_back(Mn);
    auto Fe = pNistManager->FindOrBuildElement(26);
    elements.push_back(Fe);
    auto Ni = pNistManager->FindOrBuildElement(28);
    elements.push_back(Ni);
    auto Cu = pNistManager->FindOrBuildElement(29);
    elements.push_back(Cu);
    auto Xe = pNistManager->FindOrBuildElement(54);
    elements.push_back(Xe);

    G4cout << "Avaliable elements: " << G4endl;
    for (auto pe : elements) {
        G4cout << pe->GetName() << " ";
    }
    G4cout << G4endl;

    G4double Temperature = 173.15 * kelvin;
    G4double Pressure = 1.5 * atmosphere;

    // ==== Vacuum ====
    auto Vacuum = new G4Material("Vacuum", 1.e-20 * g / cm3, 2, kStateGas);
    Vacuum->AddElement(N, 0.755);
    Vacuum->AddElement(O, 0.245);

    G4MaterialPropertiesTable *pVacuumPropertiesTable = new G4MaterialPropertiesTable();
    G4double pVacuumPhotonMomentum[] = {6.91 * eV, 6.98 * eV, 7.05 * eV};
    G4double pVacuumRefractiveIndex[] = {1.0, 1.0, 1.0};
    pVacuumPropertiesTable->AddProperty("RINDEX", pVacuumPhotonMomentum,
                                    pVacuumRefractiveIndex, 3);
    Vacuum->SetMaterialPropertiesTable(pVacuumPropertiesTable);
    materials.push_back(Vacuum);

    // ==== Air ====
    auto Air = new G4Material("Air", 0.00120479 * g / cm3, 4, kStateGas);
    Air->AddElement(C, 0.000124);
    Air->AddElement(N, 0.755268);
    Air->AddElement(O, 0.231781);
    Air->AddElement(Ar, 0.012827);

    G4MaterialPropertiesTable *pAirPropertiesTable = new G4MaterialPropertiesTable();
    G4double pAirPhotonMomentum[] = {6.91 * eV, 6.98 * eV, 7.05 * eV};
    G4double pAirRefractiveIndex[] = {1.0, 1.0, 1.0};
    pAirPropertiesTable->AddProperty("RINDEX", pAirPhotonMomentum,
                                    pAirRefractiveIndex, 3);
    Air->SetMaterialPropertiesTable(pAirPropertiesTable);
    materials.push_back(Air);

    // ==== Liquid xenon ====
    auto LXe = new G4Material("LXe", 2.862 * g / cm3, 1, kStateLiquid,
                            Temperature, Pressure);
    LXe->AddElement(Xe, 1);

    G4double pLXePhotonMomentum[] = {6.91 * eV, 6.98 * eV, 7.05 * eV};
    G4double pLXeRefractiveIndex[] = {1.63, 1.61, 1.58};
    G4double pLXeScintillation[] = {0.1, 1.0, 0.1};
    G4double pLXeAbsorbtionLength[] = {1000. * cm, 1000. * cm, 1000. * cm};
    G4double pLXeScatteringLength[] = {30. * cm, 30. * cm, 30. * cm};
    G4MaterialPropertiesTable *pLXePropertiesTable = new G4MaterialPropertiesTable();
    pLXePropertiesTable->AddProperty("RINDEX", pLXePhotonMomentum,
                                    pLXeRefractiveIndex, 3);
    pLXePropertiesTable->AddProperty("ABSLENGTH", pLXePhotonMomentum,
                                    pLXeAbsorbtionLength, 3);
    pLXePropertiesTable->AddProperty("RAYLEIGH", pLXePhotonMomentum,
                                    pLXeScatteringLength, 3);
    pLXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 3000. / keV);
    pLXePropertiesTable->AddConstProperty("RESOLUTIONSCALE", 1.0);
#if G4VERSION_NUMBER>=1100
    pLXePropertiesTable->AddProperty("SCINTILLATIONCOMPONENT1", pLXePhotonMomentum,
                                    pLXeScintillation, 3);
    pLXePropertiesTable->AddProperty("SCINTILLATIONCOMPONENT2", pLXePhotonMomentum,
                                    pLXeScintillation, 3);
    pLXePropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 3. * ns);
    pLXePropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 27. * ns);
    pLXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
    pLXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD2", 0.0);
#else
    pLXePropertiesTable->AddProperty("FASTCOMPONENT", pLXePhotonMomentum,
                                    pLXeScintillation, 3);
    pLXePropertiesTable->AddProperty("SLOWCOMPONENT", pLXePhotonMomentum,
                                    pLXeScintillation, 3);
    pLXePropertiesTable->AddConstProperty("FASTTIMECONSTANT", 3. * ns);
    pLXePropertiesTable->AddConstProperty("SLOWTIMECONSTANT", 27. * ns);
    pLXePropertiesTable->AddConstProperty("YIELDRATIO", 1.0);
#endif
    SetProperties(parameters, "LXe_", pLXePropertiesTable);
    LXe->SetMaterialPropertiesTable(pLXePropertiesTable);
    // LXe->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

    materials.push_back(LXe);

    // ==== Gaseous xenon ====
    auto GXe = new G4Material("GXe", 0.005887 * g / cm3, 1, kStateGas,
                            Temperature, Pressure);
    GXe->AddElement(Xe, 1);

    G4double pGXePhotonMomentum[] = {6.91 * eV, 6.98 * eV, 7.05 * eV};
    G4double pGXeRefractiveIndex[] = {1.00, 1.00, 1.00};
    G4double pGXeScintillation[] = {0.1, 1.0, 0.1};
    G4double pGXeAbsorbtionLength[] = {100 * m, 100 * m, 100 * m};
    G4double pGXeScatteringLength[] = {100 * m, 100 * m, 100 * m};
    G4MaterialPropertiesTable *pGXePropertiesTable =
        new G4MaterialPropertiesTable();
    pGXePropertiesTable->AddProperty("RINDEX", pGXePhotonMomentum,
                                    pGXeRefractiveIndex, 3);
    pGXePropertiesTable->AddProperty("ABSLENGTH", pGXePhotonMomentum,
                                    pGXeAbsorbtionLength, 3);
    pGXePropertiesTable->AddProperty("RAYLEIGH", pGXePhotonMomentum,
                                    pGXeScatteringLength, 3);
    pGXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 3000. / keV);
    pGXePropertiesTable->AddConstProperty("RESOLUTIONSCALE", 1.0);
#if G4VERSION_NUMBER>=1100
    pGXePropertiesTable->AddProperty("SCINTILLATIONCOMPONENT1", pGXePhotonMomentum,
                                    pGXeScintillation, 3);
    pGXePropertiesTable->AddProperty("SCINTILLATIONCOMPONENT2", pGXePhotonMomentum,
                                    pGXeScintillation, 3);
    pGXePropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 3. * ns);
    pGXePropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 27. * ns);
    pGXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
    pGXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD2", 0.0);
#else
    pGXePropertiesTable->AddProperty("FASTCOMPONENT", pGXePhotonMomentum,
                                    pGXeScintillation, 3);
    pGXePropertiesTable->AddProperty("SLOWCOMPONENT", pGXePhotonMomentum,
                                    pGXeScintillation, 3);
    pGXePropertiesTable->AddConstProperty("FASTTIMECONSTANT", 3. * ns);
    pGXePropertiesTable->AddConstProperty("SLOWTIMECONSTANT", 27. * ns);
    pGXePropertiesTable->AddConstProperty("YIELDRATIO", 1.0);
#endif
    SetProperties(parameters, "GXe_", pGXePropertiesTable);
    GXe->SetMaterialPropertiesTable(pGXePropertiesTable);
    materials.push_back(GXe);

    // ==== Quartz ====
    // ref: http://www.sciner.com/Opticsland/FS.htm
    G4Material *Quartz = new G4Material("Quartz", 2.201 * g / cm3, 2, kStateSolid);
    Quartz->AddElement(Si, 1);
    Quartz->AddElement(O, 2);

    G4double pQuartzPhotonMomentum[] = {1. * eV, 6.9 * eV, 6.91 * eV, 6.98 * eV, 7.05 * eV};
    G4double pQuartzRefractiveIndex[] = {1.50, 1.50, 1.50, 1.56, 1.60};
    G4double pQuartzAbsorbtionLength[] = {30 * m, 30 * m, 30 * m, 30 * m, 30 * m};
    G4MaterialPropertiesTable *pQuartzPropertiesTable =
        new G4MaterialPropertiesTable();
    pQuartzPropertiesTable->AddProperty("RINDEX", pQuartzPhotonMomentum,
                                        pQuartzRefractiveIndex, 5);
    pQuartzPropertiesTable->AddProperty("ABSLENGTH", pQuartzPhotonMomentum,
                                        pQuartzAbsorbtionLength, 5);
    SetProperties(parameters, "Quartz_", pQuartzPropertiesTable);
    Quartz->SetMaterialPropertiesTable(pQuartzPropertiesTable);
    materials.push_back(Quartz);

    // ==== PhotoCathodeAluminium ====
    G4Material *PhotoCathodeAluminium =
      new G4Material("PhotoCathodeAluminium", 8.00 * g / cm3, 1, kStateSolid);
    PhotoCathodeAluminium->AddElement(Al, 1);

    G4double pPhotoCathodePhotonMomentum[] = {1. * eV, 6.9 * eV, 6.91 * eV, 6.98 * eV, 7.05 * eV};
    G4double pPhotoCathodeRefractiveIndex[] = {1.50, 1.50, 1.50, 1.56, 1.60};
    G4double pPhotoCathodeAbsorbtionLength[] = {1. * nm, 1 * nm, 1. * nm, 1. * nm, 1. * nm};
    G4double pPhotoCathodeReflectivity[] = {0., 0., 0., 0., 0.};
    G4double pPhotoCathodeTransmittance[] = {1e-12, 1e-12, 1e-12, 1e-12, 1e-12};
    G4double pPhotoCathodeEsfficiency[] = {1., 1., 1., 1., 1.};
    G4MaterialPropertiesTable *pPhotoCathodePropertiesTable =
        new G4MaterialPropertiesTable();
    pPhotoCathodePropertiesTable->AddProperty(
        "RINDEX", pPhotoCathodePhotonMomentum, pPhotoCathodeRefractiveIndex, 5);
    pPhotoCathodePropertiesTable->AddProperty(
        "ABSLENGTH", pPhotoCathodePhotonMomentum, pPhotoCathodeAbsorbtionLength, 5);
    pPhotoCathodePropertiesTable->AddProperty(
        "REFLECTIVITY", pPhotoCathodePhotonMomentum, pPhotoCathodeReflectivity, 5);
    pPhotoCathodePropertiesTable->AddProperty(
        "EFFICIENCY", pPhotoCathodePhotonMomentum, pPhotoCathodeEsfficiency, 5);
    pPhotoCathodePropertiesTable->AddProperty(
        "TRANSMITTANCE", pPhotoCathodePhotonMomentum, pPhotoCathodeTransmittance, 5);
    SetProperties(parameters, "PhotoCathodeAluminium_", pPhotoCathodePropertiesTable);
    PhotoCathodeAluminium->SetMaterialPropertiesTable(pPhotoCathodePropertiesTable);
    materials.push_back(PhotoCathodeAluminium);

    // ==== Copper ====
    G4Material *Copper = new G4Material("Copper", 8.92 * g / cm3, 1, kStateSolid);
    Copper->AddElement(Cu, 1);

    G4double pCopperPhotonMomentum[] = {1.91 * eV, 6.98 * eV, 7.05 * eV};
    G4double pCopperReflectivity[] = {0.15, 0.2, 0.15};
    G4double pCopperTransmittance[] = {1e-12, 1e-12, 1e-12};
    G4MaterialPropertiesTable *pCopperPropertiesTable =
        new G4MaterialPropertiesTable();
    pCopperPropertiesTable->AddProperty("REFLECTIVITY", pCopperPhotonMomentum,
                                        pCopperReflectivity, 3);
    pCopperPropertiesTable->AddProperty("TRANSMITTANCE", pCopperPhotonMomentum,
                                        pCopperTransmittance, 3);
    SetProperties(parameters, "Copper_", pCopperPropertiesTable);
    Copper->SetMaterialPropertiesTable(pCopperPropertiesTable);
    materials.push_back(Copper);

    // ==== SS304LSteel ====
    G4Material *SS304LSteel =
      new G4Material("SS304LSteel", 8.00 * g / cm3, 5, kStateSolid);
    SS304LSteel->AddElement(Fe, 0.65);
    SS304LSteel->AddElement(Cr, 0.20);
    SS304LSteel->AddElement(Ni, 0.12);
    SS304LSteel->AddElement(Mn, 0.02);
    SS304LSteel->AddElement(Si, 0.01);

    G4double pSS304LSteelPhotonMomentum[] = {1.*eV, 5.*eV};
    G4double pSS304LSteelRefractiveIndex[] = {0.5, 0.5};
    G4double pSS304LSteelTransmittance[] = {1e-12, 1e-12};

    G4MaterialPropertiesTable *pSSPropertiesTable =
        new G4MaterialPropertiesTable();
    pSSPropertiesTable->
        AddProperty("REFLECTIVITY", pSS304LSteelPhotonMomentum, pSS304LSteelRefractiveIndex, 2);
    pSSPropertiesTable->
        AddProperty("TRANSMITTANCE", pSS304LSteelPhotonMomentum, pSS304LSteelTransmittance, 2);
    SetProperties(parameters, "SS304LSteel_", pSSPropertiesTable);
    SS304LSteel->SetMaterialPropertiesTable(pSSPropertiesTable);
    materials.push_back(SS304LSteel);

    // ==== Flange ====
    G4Material *Flange =
      new G4Material("Flange", 8.00 * g / cm3, 5, kStateSolid);
    Flange->AddElement(Fe, 0.65);
    Flange->AddElement(Cr, 0.20);
    Flange->AddElement(Ni, 0.12);
    Flange->AddElement(Mn, 0.02);
    Flange->AddElement(Si, 0.01);

    G4double pFlangePhotonMomentum[] = {1.*eV, 5.*eV};
    G4double pFlangeRefractiveIndex[] = {0.5, 0.5};
    G4double pFlangeTransmittance[] = {1e-12, 1e-12};

    G4MaterialPropertiesTable *pFlangePropertiesTable =
        new G4MaterialPropertiesTable();
    pFlangePropertiesTable->
        AddProperty("REFLECTIVITY", pFlangePhotonMomentum, pFlangeRefractiveIndex, 2);
    pFlangePropertiesTable->
        AddProperty("TRANSMITTANCE", pFlangePhotonMomentum, pFlangeTransmittance, 2);
    SetProperties(parameters, "Flange_", pFlangePropertiesTable);
    Flange->SetMaterialPropertiesTable(pFlangePropertiesTable);
    materials.push_back(Flange); 

    // ==== Teflon ====
    G4Material *Teflon = new G4Material("Teflon", 2.2 * g / cm3, 2, kStateSolid);
    Teflon->AddElement(C, 0.240183);
    Teflon->AddElement(F, 0.759817);

    G4double pTeflonPhotonMomentum[] = {6.91 * eV, 6.98 * eV, 7.05 * eV};
    G4double pTeflonRefractiveIndex[] = {1.63, 1.61, 1.58};
    G4double pTeflonAbsorbtionLength[] = {0.1 * cm, 0.1 * cm, 0.1 * cm};
    G4double pTeflonReflectivity[] = {0.99, 0.99, 0.99};
    G4double pTeflonSpecularLobe[] = {0.01, 0.01, 0.01};
    G4double pTeflonSpecularSpike[] = {0.01, 0.01, 0.01};
    G4double pTeflonBackscatter[] = {0.01, 0.01, 0.01};
    G4double pTeflonTransmittance[] = {1e-12, 1e-12, 1e-12};
    G4MaterialPropertiesTable *pTeflonPropertiesTable =
        new G4MaterialPropertiesTable();
    pTeflonPropertiesTable->AddProperty("RINDEX", pTeflonPhotonMomentum,
                                        pTeflonRefractiveIndex, 3);
    pTeflonPropertiesTable->AddProperty("ABSLENGTH", pTeflonPhotonMomentum,
                                        pTeflonAbsorbtionLength, 3);
    pTeflonPropertiesTable->AddProperty("REFLECTIVITY", pTeflonPhotonMomentum,
                                        pTeflonReflectivity, 3);
    pTeflonPropertiesTable->AddProperty("SPECULARLOBECONSTANT",
                                        pTeflonPhotonMomentum,
                                        pTeflonSpecularLobe, 3);
    pTeflonPropertiesTable->AddProperty("SPECULARSPIKECONSTANT",
                                        pTeflonPhotonMomentum,
                                        pTeflonSpecularSpike, 3);
    pTeflonPropertiesTable->AddProperty("BACKSCATTERCONSTANT",
                                        pTeflonPhotonMomentum,
                                        pTeflonBackscatter, 3);
    pTeflonPropertiesTable->AddProperty("TRANSMITTANCE", pTeflonPhotonMomentum,
                                        pTeflonTransmittance, 3);
    SetProperties(parameters, "Teflon_", pTeflonPropertiesTable);
    Teflon->SetMaterialPropertiesTable(pTeflonPropertiesTable);
    materials.push_back(Teflon);

    // ==== Cirlex ====
    G4Material *Cirlex = new G4Material("Cirlex", 1.43 * g / cm3, 4, kStateSolid);
    Cirlex->AddElement(C, 22);
    Cirlex->AddElement(H, 10);
    Cirlex->AddElement(N, 2);
    Cirlex->AddElement(O, 5);

    G4double pCirlexPhotonMomentum[] = {6.91*eV, 6.98*eV, 7.05*eV};
    G4double pCirlexReflectivity[] = {0.5, 0.5, 0.5};
    G4double pCirlexTransmittance[] = {1e-12, 1e-12, 1e-12};
    G4MaterialPropertiesTable *pCirlexPropertiesTable =
        new G4MaterialPropertiesTable();
    pCirlexPropertiesTable->AddProperty("REFLECTIVITY", pCirlexPhotonMomentum,
                                        pCirlexReflectivity, 3);
    pCirlexPropertiesTable->AddProperty("TRANSMITTANCE", pCirlexPhotonMomentum,
                                        pCirlexTransmittance, 3);
    SetProperties(parameters, "Cirlex_", pCirlexPropertiesTable);
    Cirlex->SetMaterialPropertiesTable(pCirlexPropertiesTable);
    materials.push_back(Cirlex);

    // ==== Polyethylene ====
    // https://www.sciencedirect.com/science/article/pii/S0969806X18311289
    // Fig3. 10% Boron -> 1.03g/cm^3
    G4Material *Polyethylene = new G4Material("Polyethylene", 1.03 * g / cm3, 3, kStateSolid);
    Polyethylene->AddElementByMassFraction(C, 0.779);
    Polyethylene->AddElementByMassFraction(H, 0.121);
    Polyethylene->AddElementByMassFraction(B, 0.1);

    // ==== Nist material ====
    auto air = pNistManager->FindOrBuildMaterial("G4_AIR");
    materials.push_back(air);

    auto polyEthylene = pNistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
    materials.push_back(polyEthylene);

    // auto teflon = pNistManager->FindOrBuildMaterial("G4_TEFLON");
    // materials.push_back(teflon);

    // auto copper = pNistManager->FindOrBuildMaterial("G4_Cu");
    // materials.push_back(copper);

    auto lead = pNistManager->FindOrBuildMaterial("G4_Pb");
    materials.push_back(lead);

    // auto liquidXenon = pNistManager->FindOrBuildMaterial("G4_lXe");
    // materials.push_back(liquidXenon);

    // auto stainlessSteel =
    //     pNistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    // materials.push_back(stainlessSteel);

    // ==== concrete ====
    auto concrete = new G4Material("concrete", 2.3 * g / cm3, 4, kStateSolid);
    // 混凝土的成分
    concrete->AddElement(G4Element::GetElement("O"), 0.52);   // 氧的质量分数
    concrete->AddElement(G4Element::GetElement("Si"), 0.25);  // 硅的质量分数
    concrete->AddElement(G4Element::GetElement("Ca"), 0.17);  // 钙的质量分数
    concrete->AddElement(G4Element::GetElement("Al"), 0.06);  // 铝的质量分数
    materials.push_back(concrete);
    //数字是整数则是数量比，是浮点数且加起来为1则为质量比，浮点数但相加不为一则自动缩放，但可能有bug
    // ==== water ====
    auto water = pNistManager->FindOrBuildMaterial("G4_WATER");
    materials.push_back(water);

    // ==== PEEK ====
    auto PEEK = new G4Material("PEEK", 1.32 * g / cm3, 3, kStateSolid);
    PEEK->AddElement(C, 19);
    PEEK->AddElement(H, 12);
    PEEK->AddElement(O, 3);
    materials.push_back(PEEK);

    // ==== PMTbasePEEK ====
    auto PMTbasePEEK = new G4Material("PMTbasePEEK", 1.32 * g / cm3, 3, kStateSolid);
    PMTbasePEEK->AddElement(C, 19);
    PMTbasePEEK->AddElement(H, 12);
    PMTbasePEEK->AddElement(O, 3);
    materials.push_back(PMTbasePEEK);


    // EJ200塑料闪烁体（PS)
    auto PS = new G4Material("PS", 1.023 * g / cm3, 2, kStateSolid);
    PS->AddElement(C, 0.915);
    PS->AddElement(H, 0.085);
    materials.push_back(PS);
    // 质量分数来源：https://eljentechnology.com/products/plastic-scintillators/ej-200-ej-204-ej-208-ej-212

    G4cout << "Available materials: " << G4endl;
    for (auto pm : materials) {
        G4cout << pm->GetName() << " ";
    }
    G4cout << G4endl;
}

// Set additional properties for material based on configuration file
// User should define parameter(s) in .xml or .json's 'material' element
// 
// Example:
// <material name="Material">
//   <parameter name="LXe_SCINTILLATIONYIELD" value="3000000."/>
//   <parameter name="Teflon_RINDEX" value="6.91 1.63 6.98 1.61 7.05 1.58"/>
// </material>
void Material::SetProperties(
    const BambooParameters &parameters,
    const std::string &material,
    G4MaterialPropertiesTable *table
) {
    bool changed_flag = false;
    for (std::string property : vector_properties) {
        auto pstring = parameters.getParameter(material + property);
        if (pstring != "") {
            // G4cout << pstring << G4endl;
            SetVectorProperty(table, property, pstring);
            changed_flag = true;
        }
    }

    for (std::string property : const_properties) {
        auto pstring = parameters.getParameter(material + property);
        if (pstring != "") {
            // G4cout << pstring << G4endl;
            SetConstProperty(table, property, pstring);
            changed_flag = true;
        }
    }

    if (changed_flag) {
        G4cout << "The MPT for the " + material + " is now: " << G4endl;
        table->DumpTable();
        G4cout << "............." << G4endl;
    }
}

// Define property as function of energy based on configuration file
// User should set parameters for 'material' element
// Energy in configuration file should be in eV
// Energy and property value should be written alternately
void Material::SetVectorProperty(
    G4MaterialPropertiesTable *table,
    const std::string &property,
    G4String newValue
) {
    G4MaterialPropertyVector* mpv = new G4MaterialPropertyVector();
    std::istringstream instring(newValue);
    while(instring)
    {
        G4String tmp;
        instring >> tmp;
        if ( tmp == "" ) { break; }
        // energy in eV!
        auto en = G4UIcommand::ConvertToDouble(tmp) / 1e6;
        instring >> tmp;
        auto val = G4UIcommand::ConvertToDouble(tmp);
        mpv->InsertValues(en, val);
    }

    table->AddProperty(property.c_str(), mpv);
}

// Define constant property based on configuration file
void Material::SetConstProperty(
    G4MaterialPropertiesTable *table,
    const std::string &property,
    G4String newValue
) {
    std::istringstream instring(newValue);
    G4String tmp;
    instring >> tmp;
    auto val = G4UIcommand::ConvertToDouble(tmp);

    table->AddConstProperty(property.c_str(), val);
}

// Evaluate optical parameters based on configuration file
// The optional properties are in private member variables of Material
// 
// Example:
// <parameter name="base_type" value="dielectric_metal"/>
// <parameter name="casing_finish" value="polished"/>
// <parameter name="casing_type" value="dielectric_metal"/>
// <parameter name="window_finish" value="polished"/>
// <parameter name="photocathode_finish" value="polished"/>
G4OpticalSurface* Material::EvalOpParameters(
    const BambooParameters &parameters,
    G4VPhysicalVolume* mainPV,
    const G4String &prefix
) {
    auto material = parameters.getParameter(prefix + "material");
    if (material == "") {
        material = mainPV->GetLogicalVolume()->GetMaterial()->GetName();
    }
    auto model = parameters.getParameter(prefix + "model");
    auto finish = parameters.getParameter(prefix + "finish");
    auto type = parameters.getParameter(prefix + "type");
    G4double alpha = 0.1;  // default value
    if (parameters.getParameters().find(prefix + "alpha") != parameters.getParameters().end()) {
        alpha = parameters.evaluateParameter(prefix + "alpha");
    }
    auto medium = G4Material::GetMaterial(material);
    auto old_table = medium->GetMaterialPropertiesTable();
    // need to copy-construct a new table!
    G4MaterialPropertiesTable *new_table = new G4MaterialPropertiesTable(*old_table);
    // update properties of surface
    SetProperties(parameters, prefix, new_table);

    G4OpticalSurface* opSurface = new G4OpticalSurface(
        material + "_Surface",
        Models.at(model),
        Finishes.at(finish),
        Types.at(type),
        alpha
    );
    opSurface->SetMaterialPropertiesTable(new_table);
    return opSurface;
}

// Add LogicalSkinSurface for mono-material block
// No need to input parent volume
void Material::AddLogicalSkinSurface(
    G4LogicalVolume* mainLV,
    G4OpticalSurface* opSurface,
    const BambooParameters &parameters,
    const G4String &prefix
) {
    auto material = parameters.getParameter(prefix + "material");

    G4Material* medium;
    if (material != "") {
        medium = G4Material::GetMaterial(material);
    } else {
        medium = mainLV->GetMaterial();
    }

    new G4LogicalSkinSurface(medium->GetName() + "_Skin", mainLV, opSurface);
    G4cout << "Set G4LogicalSkinSurface to "\
    "" + mainLV->GetName() + " with " + opSurface->GetName() << G4endl;
}

// Add LogicalBorderSurface for border of PhysicalVolumes
// If border is between volume and its sister volume
// (which means not between volume and its mother volume),
// search(deduce) sister volume by name based on input mother volume's PhysicalVolume
void Material::AddLogicalBorderSurface(
    G4VPhysicalVolume* mainPV,
    G4VPhysicalVolume* containerPV,
    G4OpticalSurface* opSurface,
    const BambooParameters &parameters,
    const G4String &prefix
) {
    G4String partner = parameters.getParameter(prefix + "partner");
    auto order = parameters.evaluateParameter(prefix + "order");

    BorderPVs borderpvs;
    if (partner == "parent") {
        if (order == 0) {
            borderpvs.vol1 = mainPV;
            borderpvs.vol2 = containerPV;
        } else {
            borderpvs.vol1 = containerPV;
            borderpvs.vol2 = mainPV;
        }
    } else {
        auto containerLV = containerPV->GetLogicalVolume();
        auto nDaughters = containerLV->GetNoDaughters();
        G4VPhysicalVolume* sisterPV;
        int id;
        for (id = 0; id < nDaughters; ++id) {
            sisterPV = containerLV->GetDaughter(id);
            if (sisterPV->GetName() == partner) {
                break;
            }
        }
        if (id == nDaughters) {
            throw std::invalid_argument("sister physical volume " + partner + " not defined");
        }

        if (order == 0) {
            borderpvs.vol1 = mainPV;
            borderpvs.vol2 = sisterPV;
        } else {
            borderpvs.vol1 = sisterPV;
            borderpvs.vol2 = mainPV;
        }
    }

    new G4LogicalBorderSurface(
        borderpvs.vol1->GetName() + "_to_" + borderpvs.vol2->GetName() + "_Border",
        borderpvs.vol1, borderpvs.vol2, opSurface
    );
    G4cout << "Set G4LogicalBorderSurface from "\
    "" + borderpvs.vol1->GetName() + " to "\
    "" + borderpvs.vol2->GetName() + " with " + opSurface->GetName() << G4endl;
}

// Define logical surface based on parameter(s) of LV and border of PVs
// If defined partner(daughter volume or sister volume) in component's parameter, 
// add LogicalBorderSurface, otherwise, add LogicalSkinSurface.
void Material::AddLogicalSurface(
    G4VPhysicalVolume* mainPV,
    G4VPhysicalVolume* containerPV,
    const BambooParameters &parameters,
    const G4String &prefix
) {
    auto opSurface = EvalOpParameters(parameters, mainPV, prefix);

    auto partner = parameters.getParameter(prefix + "partner");
    if (partner == "") {
        // if not sister volume name provided,
        // just assign LogicalSkinSurface to LV
        AddLogicalSkinSurface(
            mainPV->GetLogicalVolume(),
            opSurface,
            parameters,
            prefix
        );
    } else {
        AddLogicalBorderSurface(
            mainPV,
            containerPV,
            opSurface,
            parameters,
            prefix
        );
    }
}

