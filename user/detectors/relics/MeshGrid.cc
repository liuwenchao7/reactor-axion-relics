//hexagon meshgrid
#include "MeshGrid.hh"
#include "PandaXSensitiveDetector.hh"

#include <G4Version.hh>
#include <G4TwoVector.hh>
#include <G4Material.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4LogicalVolume.hh>
#include <G4UnionSolid.hh>
#include <G4ExtrudedSolid.hh>
#include <G4Polyhedra.hh>
#include <G4GenericPolycone.hh>
#include <G4SubtractionSolid.hh>
#if G4VERSION_NUMBER>=1100
  #include <G4MultiUnion.hh>
#endif

DetectorRegister<MeshGrid, std::string, BambooParameters> MeshGrid::reg("MeshGrid");

MeshGrid::MeshGrid (const std::string &n, const BambooParameters &pars)
  : Monoblock(n, pars) {
    G4cout << "create detector MeshGrid..." << G4endl;
}

bool MeshGrid::constructMainLV(const BambooParameters &) {
    using namespace CLHEP;
    auto soliname = parameters.getParameter("soliname");
    auto logivol = parameters.getParameter("logivol");
    std::string source = parameters.getParameter("source");
    auto radius = parameters.evaluateParameter("radius");
    auto shape = parameters.getParameter("shape");
    auto wire_diameter = parameters.evaluateParameter("wire_diameter");
    auto wire_pitch = parameters.evaluateParameter("wire_pitch");//本义网格间距，本代码作为六边形边长使用
    // default value of sides is 0
    int sides = parameters.getParameter<int>("sides");
    auto phi = parameters.evaluateParameter("phi");
    int fold = parameters.evaluateParameter("fold");
    auto ring_width = parameters.evaluateParameter("ring_width");
    auto ring_height = parameters.evaluateParameter("ring_height");
    auto material = parameters.getParameter("material");
    auto optical_sister = parameters.getParameter("optical_sister");
    if (radius == 0) {
        radius = 1 * m;
    }
    if (shape == "") {
        shape = "Box";
    }
    if (wire_diameter == 0.) {
        wire_diameter = 0.2 * mm;
    }
    if (wire_pitch == 0.) {
        wire_pitch = 2.89 * mm;
    }
    if (fold == 0) {
        fold = 2;
    } else if (fold > 2) {
        throw std::invalid_argument("Too large fold, which is " + fold);
    }


    //以上为没有传参的情况下的默认设置

    G4double inner;
    G4double outer;
    if (sides == 0) {
        inner = radius;
        outer = radius + ring_width;
    } else {
        inner = radius / std::cos(M_PI / sides);
        outer = (radius + ring_width) / std::cos(M_PI / sides);
    }
    int n_wires_x = floor((radius / mm) / (sqrt(3) * wire_pitch / mm));
    int n_wires_z = floor((radius / mm) / (1.5 * wire_pitch / mm));

#if G4VERSION_NUMBER>=1100
    auto electrode = new G4MultiUnion("electrode");
#else
    auto small_box = new G4Box(
        "small_box", wire_diameter / 100, wire_diameter / 100, wire_diameter / 100);//基点
    auto electrode = new G4UnionSolid(soliname, small_box, small_box);
#endif
    G4VSolid* wire;//网格基本单元wire
    if (shape == "Tubs") {
        wire = new G4Tubs(soliname, 0, wire_diameter / 2, wire_pitch , 0, 2 * M_PI);
    } else if (shape == "Box") {
        wire = new G4Box(soliname, wire_diameter / 2, wire_diameter / 2,  wire_pitch/2 );
    } else {
        throw std::invalid_argument("Wire shape " + shape + " not recognized");
    }

    for (int jl = -n_wires_z - 1; jl <= n_wires_z + 1; ++jl) {
        for (int jp = -n_wires_x - 1; jp <= n_wires_x + 1; ++jp){
            G4ThreeVector position1 = G4ThreeVector((jp + abs(jl % 2)/2.0) * sqrt(3) * wire_pitch, 0., 1.5 * jl * wire_pitch);
            G4Transform3D transform1 = G4Transform3D(G4RotationMatrix(), position1);//竖直部分
#if     G4VERSION_NUMBER>=1100//版本号！！适应不同版本
            electrode->AddNode(wire, transform1);
#else
            electrode = new G4UnionSolid(soliname, electrode, wire, transform1);
#endif

            G4ThreeVector position2 = G4ThreeVector((jp + abs(jl % 2)/2.0) * sqrt(3) * wire_pitch + sqrt(3)/4.0 * wire_pitch, 0., 1.5 * jl * wire_pitch - 3/4.0 * wire_pitch );
            G4RotationMatrix *Rotation2 = new G4RotationMatrix();
            Rotation2->rotateY(-60 * deg);
            G4Transform3D transform2 = G4Transform3D(*Rotation2, position2);
#if     G4VERSION_NUMBER>=1100//版本号！！适应不同版本
            electrode->AddNode(wire, transform2);
#else
            electrode = new G4UnionSolid(soliname, electrode, wire, transform2);
#endif

            G4ThreeVector position3 = G4ThreeVector((jp + abs(jl % 2)/2.0) * sqrt(3) * wire_pitch + sqrt(3)/4.0 * wire_pitch, 0., 1.5 * jl * wire_pitch + 3/4.0 * wire_pitch );
            G4RotationMatrix *Rotation3 = new G4RotationMatrix();
            Rotation3->rotateY(60 * deg);
            G4Transform3D transform3 = G4Transform3D(*Rotation3, position3);
#if     G4VERSION_NUMBER>=1100//版本号！！适应不同版本
            electrode->AddNode(wire, transform3);
#else
            electrode = new G4UnionSolid(soliname, electrode, wire, transform3);
#endif
             
        }
       
    }
#if G4VERSION_NUMBER>=1100
        electrode->Voxelize();
#endif

    auto origin = new G4Box(
        "origin", wire_diameter / 10, wire_diameter / 10, wire_diameter / 10);//origin为基点

    G4RotationMatrix *mRot = new G4RotationMatrix();
    mRot->rotateX(90 * deg);//旋转后网格法向为z方向
    mRot->rotateZ(phi * deg);//phi
    G4Transform3D transform = G4Transform3D(*mRot, G4ThreeVector());
    auto meshgrid = new G4UnionSolid(soliname, origin, electrode, transform);//并集，transform 定义了 electrode 相对于 origin 的位置和方向。

    if (ring_width != 0) {//ring??
        G4VSolid* ring_solid;
        if (sides == 0) {
            ring_solid = new G4Tubs(soliname, inner, 1.1 * outer, ring_height / 2, 0, 2 * M_PI);
        } else {
            G4double ring_r[4] = {inner, 1.1 * outer, 1.1 * outer, inner};
            G4double ring_z[4] = {-ring_height / 2, -ring_height / 2, ring_height / 2, ring_height / 2};
            ring_solid = new G4Polyhedra(soliname, 0, 2 * M_PI, sides, 4, ring_r, ring_z);
        }
        meshgrid = new G4UnionSolid(soliname, ring_solid, meshgrid);
    }

    G4double largest_height = std::max(ring_height / mm, wire_diameter / mm) * mm;
    G4double large_cover = 2 * outer;
    G4VSolid* veto_solid;
    if (sides == 0) {
        veto_solid = new G4Tubs(soliname, outer, large_cover, largest_height, 0, 2 * M_PI);
    } else {
        G4TwoVector offsetA(0, 0), offsetB(0, 0);//二维向量
        G4double scaleA = 1, scaleB = 1;
        veto_solid = new G4Tubs(soliname, 0, large_cover, largest_height, 0, 2 * M_PI);
        auto polygon = GetPolygon(sides, outer);
        auto veto_veto_solid = new G4ExtrudedSolid(
            soliname, polygon, 2 * largest_height, offsetA, scaleA, offsetB, scaleB);
        veto_solid = new G4SubtractionSolid(soliname, veto_solid, veto_veto_solid);
    }

    auto solid = new G4SubtractionSolid(soliname, meshgrid, veto_solid);
    auto medium = G4Material::GetMaterial(material);
    mainLV =
        new G4LogicalVolume(solid, medium, logivol, 0, 0, 0);
    return true;
}

std::vector<G4TwoVector> MeshGrid::GetPolygon(
    const G4int &sides,
    const G4double &radius
) {
    std::vector<G4TwoVector> polygon(sides);
    G4double angle = 2 * M_PI / sides;
    for (G4int il = 0; il < sides; ++il)
    {
        G4double phii = il * angle;
        G4double cosphii = std::cos(phii);
        G4double sinphii = std::sin(phii);
        polygon[il].set(radius * cosphii, radius * sinphii);
    }
    return polygon;
}