#include "PMTsR8520.hh"
#include "Material.hh"

#include <G4Material.hh>
#include <G4ThreeVector.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4VisAttributes.hh>
#include <G4Box.hh>
#include <G4Polyhedra.hh>
#include <G4UnionSolid.hh>
#include <G4Transform3D.hh>
#include <G4SystemOfUnits.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalBorderSurface.hh>

PMTsR8520::PMTsR8520 (const std::string &, const BambooParameters &pars)
    : parameters{pars} {

    G4cout << "create detector PMTsR8520..." << G4endl;
    WindowWidth = 25.00 * mm;
    WindowThickness = 1.50 * mm;

    PhotoCathodeWidth = 22.00 * mm;
    PhotoCathodeThickness = 0.50 * mm;

    CasingHeight = 27.00 * mm;
    CasingWidth = 25.40 * mm;
    CasingThickness = 0.50 * mm;

    BaseThickness = 1.50 * mm;
    PMTToPMTBase = 3.00 * mm;
}

PMTLVs PMTsR8520::construct (const BambooParameters &, BambooDetector *parent) {
    // add construction code here
    using namespace CLHEP;
    G4Material *Vacuum = G4Material::GetMaterial("Vacuum");
    G4Material *Quartz = G4Material::GetMaterial("Quartz");
    G4Material *Cirlex = G4Material::GetMaterial("Cirlex");
    G4Material *PhotoCathodeAluminium =
        G4Material::GetMaterial("PhotoCathodeAluminium");
    G4Material *SS304LSteel = G4Material::GetMaterial("SS304LSteel");

    bool check_overlap = parameters.getParameter<bool>("check_overlap");
    auto parent_medium = parent->getMainPV()->GetLogicalVolume()->GetMaterial();

    Material* pMaterial = new Material();

    const G4double HalfZ =
        0.5 * (WindowThickness + CasingHeight + PMTToPMTBase + BaseThickness);

    // Base
    const G4double BaseHalfX = 0.5 * CasingWidth;
    const G4double BaseHalfZ = 0.5 * BaseThickness;
    const G4double BaseOffset = HalfZ - BaseHalfZ;
    G4Box *BaseBox = new G4Box(
        "BaseBox",
        BaseHalfX,
        BaseHalfX,
        BaseHalfZ
    );

    BaseLogicalVolume = new G4LogicalVolume(
        BaseBox, Cirlex, "BaseLogicalVolume", 0, 0, 0);

    // Window
    const G4double WindowHalfX = 0.5 * WindowWidth;
    const G4double WindowHalfZ = 0.5 * WindowThickness;
    const G4double WindowOffset = -HalfZ + WindowHalfZ;
    G4Box *WindowBox = new G4Box(
        "WindowBox",
        WindowHalfX,
        WindowHalfX,
        WindowHalfZ
    );

    WindowLogicalVolume = new G4LogicalVolume(
        WindowBox, Quartz, "WindowLogicalVolume", 0, 0, 0);

    // Casing
    const G4double CasingHalfX = 0.5 * CasingWidth;
    const G4double CasingHalfZ = 0.5 * CasingHeight;
    const G4double CasingOffset = -HalfZ + WindowThickness + CasingHalfZ;

    G4Box *CasingBox = new G4Box(
        "CasingBox",
        CasingHalfX,
        CasingHalfX,
        CasingHalfZ
    );

    CasingLogicalVolume = new G4LogicalVolume(
        CasingBox, SS304LSteel, "CasingLogicalVolume", 0, 0, 0);

    // Interior
    const G4double InteriorHalfX = CasingHalfX - CasingThickness;
    const G4double InteriorHalfZ = CasingHalfZ - CasingThickness;

    G4Box *InteriorBox = new G4Box(
        "InteriorBox", InteriorHalfX,
        InteriorHalfX, InteriorHalfZ);

    InteriorLogicalVolume = new G4LogicalVolume(
        InteriorBox, Vacuum, "InteriorLogicalVolume", 0, 0, 0);

    // Photocathode
    const G4double PhotoCathodeHalfX = 0.5 * PhotoCathodeWidth;
    const G4double PhotoCathodeHalfZ = 0.5 * PhotoCathodeThickness;
    const G4double PhotoCathodeOffsetZ =
        PhotoCathodeHalfZ - CasingHalfZ;

    G4Box *PhotoCathodeBox =
        new G4Box(
            "PhotoCathodeBox", PhotoCathodeHalfX,
            PhotoCathodeHalfX, PhotoCathodeHalfZ);

    PhotoCathodeLogicalVolume =
        new G4LogicalVolume(PhotoCathodeBox, PhotoCathodeAluminium,
                            "PhotoCathodeLogicalVolume", 0, 0, 0);

    // Combine
    G4UnionSolid* CombineBox;
    G4Transform3D transform;
    auto origin = new G4Box("origin", 1. * mm, 1. * mm, 1. * mm);
    auto mRot = G4RotationMatrix();
    G4ThreeVector BaseTrans(0, 0, BaseOffset);
    transform = G4Transform3D(mRot, BaseTrans);
    CombineBox = new G4UnionSolid("CombineBox", origin, BaseBox, transform);
    G4ThreeVector CasingTrans(0, 0, CasingOffset);
    transform = G4Transform3D(mRot, CasingTrans);
    CombineBox = new G4UnionSolid("CombineBox", CombineBox, CasingBox, transform);
    G4ThreeVector WindowTrans(0, 0, WindowOffset);
    transform = G4Transform3D(mRot, WindowTrans);
    CombineBox = new G4UnionSolid("CombineBox", CombineBox, WindowBox, transform);

    CombineLogicalVolume =
        new G4LogicalVolume(CombineBox, parent_medium, "CombineLogicalVolume", 0, 0, 0);

    // Set optical surfaces
    BasePhysicalVolume = new G4PVPlacement(
        0, G4ThreeVector(0., 0., BaseOffset), BaseLogicalVolume,
        "Base", CombineLogicalVolume, false, 0, check_overlap);

    WindowPhysicalVolume = new G4PVPlacement(
        0, G4ThreeVector(0., 0., WindowOffset), WindowLogicalVolume,
        "Window", CombineLogicalVolume, false, 0, check_overlap);

    CasingPhysicalVolume = new G4PVPlacement(
        0, G4ThreeVector(0., 0., CasingOffset), CasingLogicalVolume,
        "Casing", CombineLogicalVolume, false, 0, check_overlap);

    InteriorPhysicalVolume = new G4PVPlacement(
        0, G4ThreeVector(0., 0., 0.), InteriorLogicalVolume,
        "Interior", CasingLogicalVolume, false, 0, check_overlap);

    PhotoCathodePhysicalVolume = new G4PVPlacement(
        0, G4ThreeVector(0., 0., PhotoCathodeOffsetZ),
        PhotoCathodeLogicalVolume, "PhotoCathode",
        CasingLogicalVolume, false, 0, check_overlap);

    // Set optical surfaces
    pMaterial->AddLogicalSurface(BasePhysicalVolume, parent->getMainPV(), parameters, "base_");
    pMaterial->AddLogicalSurface(WindowPhysicalVolume, parent->getMainPV(), parameters, "window_");
    pMaterial->AddLogicalSurface(CasingPhysicalVolume, parent->getMainPV(), parameters, "casing_");
    pMaterial->AddLogicalSurface(PhotoCathodePhysicalVolume, parent->getMainPV(), parameters, "photocathode_");

    // Attributes
    G4Colour hBaseColor(0.8, 0.4, 0.1, 0.2);
    G4VisAttributes *BaseVisAtt = new G4VisAttributes(hBaseColor);
    BaseVisAtt->SetVisibility(true);
    BaseLogicalVolume->SetVisAttributes(BaseVisAtt);

    G4Colour hWindowColor(0., 1., 1., 0.5);
    G4VisAttributes *WindowVisAtt = new G4VisAttributes(hWindowColor);
    WindowVisAtt->SetVisibility(true);
    WindowLogicalVolume->SetVisAttributes(WindowVisAtt);

    G4Colour hPhotocathodeColor(1., 0.5, 0.1, 0.9);
    G4VisAttributes *PhotocathodeVisAtt =
        new G4VisAttributes(hPhotocathodeColor);
    PhotocathodeVisAtt->SetVisibility(true);
    PhotoCathodeLogicalVolume->SetVisAttributes(PhotocathodeVisAtt);

    G4Colour hColor(1., 1., 1., 0.7);
    G4VisAttributes *VisAtt = new G4VisAttributes(hColor);
    VisAtt->SetVisibility(true);
    CombineLogicalVolume->SetVisAttributes(VisAtt);

    PMTLVs pmtlvs;
    pmtlvs.CombineLV = CombineLogicalVolume;
    pmtlvs.WindowLV = WindowLogicalVolume;
    pmtlvs.PhotoCathodeLV = PhotoCathodeLogicalVolume;
    pmtlvs.CasingLV = CasingLogicalVolume;

    return pmtlvs;
}

