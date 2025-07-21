#include <STEPControl_Reader.hxx> 
#include <IMeshTools_Parameters.hxx>
#include <BRep_Tool.hxx>
#include <BRepLib.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <iostream>
#include <fstream>
void readMultiPiece(const float deflection_err, TopoDS_Shape& shape,
                      std::vector<std::vector<float>>& points,
                      std::vector<std::vector<uint32_t>>& indices) {
  // Triangulation parameters.
  IMeshTools_Parameters aMeshParams;
  aMeshParams.Deflection = deflection_err;
  aMeshParams.Angle = 0.5;
  aMeshParams.Relative = Standard_False;
  aMeshParams.InParallel = Standard_True;
  aMeshParams.MinSize = Precision::Confusion();
  aMeshParams.InternalVerticesMode = Standard_True;
  aMeshParams.ControlSurfaceDeflection = Standard_True;
  aMeshParams.DeflectionInterior = deflection_err;

  BRepMesh_IncrementalMesh(shape, aMeshParams);

  TopExp_Explorer partExplorer;

  for (partExplorer.Init(shape, TopAbs_SOLID); partExplorer.More();
       partExplorer.Next()) {
    TopoDS_Shape part = partExplorer.Current();

    uint32_t trigIdxOffset = 0;
    std::vector<float> pointsCur;
    std::vector<uint32_t> indicesCur;

    for (TopExp_Explorer aExpFace(part, TopAbs_FACE); aExpFace.More();
         aExpFace.Next()) {
      TopoDS_Face aFace = TopoDS::Face(aExpFace.Current());
      TopAbs_Orientation faceOrientation = aFace.Orientation();

      TopLoc_Location aLocation;
      Handle(Poly_Triangulation) aTr =
          BRep_Tool::Triangulation(aFace, aLocation);

      if (!aTr.IsNull()) {
        const Poly_Array1OfTriangle& triangles = aTr->Triangles();
        TColgp_Array1OfPnt aPoints(1, aTr->NbNodes());
        for (Standard_Integer i = 1; i <= aTr->NbNodes(); i++)
          aPoints(i) = aTr->Node(i).Transformed(aLocation);

        // insert points to array.
        for (int i = 1; i <= aPoints.Size(); i++) {
          pointsCur.push_back(aPoints[i].X());
          pointsCur.push_back(aPoints[i].Y());
          pointsCur.push_back(aPoints[i].Z());
        }

        Standard_Integer nF = aTr->NbTriangles();
        Standard_Integer nt, n1, n2, n3;
        for (nt = 1; nt <= nF; nt++) {
          triangles(nt).Get(n1, n2, n3);
          if (faceOrientation != TopAbs_Orientation::TopAbs_FORWARD) {
            std::swap(n1, n3);
          }
          // c++ indices start with 0, occt indices start with 1.
          indicesCur.push_back(n1 + trigIdxOffset - 1);
          indicesCur.push_back(n2 + trigIdxOffset - 1);
          indicesCur.push_back(n3 + trigIdxOffset - 1);
        }
        trigIdxOffset += aPoints.Size();
      }
    }

    points.push_back(pointsCur);
    indices.push_back(indicesCur);
  }
}

int main(int argc, char*argv[]){
  STEPControl_Reader    reader;
  const std::string     fileName = "F:/meshes/step/roundCube.step";
  IFSelect_ReturnStatus ret = reader.ReadFile(fileName.c_str());
  if (ret != IFSelect_RetDone) {  // read fail.
    std::cout << "OCCT step read Fail" << std::endl;
    return false;
  }
  Standard_Integer NbRoots = reader.NbRootsForTransfer();
  Standard_Integer NbTrans = reader.TransferRoots();
  TopoDS_Shape shape = reader.OneShape();
  std::vector<std::vector<float>> points ;
  std::vector<std::vector<uint32_t>> indices;
  float deflectionError = 0.01;
  readMultiPiece(deflectionError, shape, points, indices);
  std::cout << points.size() << "\n";
  std::cout << indices.size() << "\n";


  for (size_t i = 0; i < points.size(); i++) {
    std::string outFile = fileName + "_" + std::to_string(i) + ".obj";
    std::ofstream out(outFile);
    const std::vector<float>& v = points[i];
    for (size_t j = 0; j < v.size(); j += 3) {
      out << "v ";
      out << v[j] << " " << v[j + 1] << " " << v[j + 2] << "\n";
    }
    const std::vector<uint32_t>& t = indices[i];
    for (size_t j = 0; j < t.size(); j += 3) {
      out << "f ";
      out << (t[j] + 1) << " " << (t[j + 1] + 1) << " " << (t[j + 2] + 1)
          << "\n";
    }
    out.close();
  }
	return 0;
}