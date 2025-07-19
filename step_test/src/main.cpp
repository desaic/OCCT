#include <STEPControl_Reader.hxx> 

int main(int argc, char*argv[]){
  STEPControl_Reader    reader;
  const std::string     fileName = "F:/meshes/step/roundCube.step";
  IFSelect_ReturnStatus ret = reader.ReadFile(fileName.c_str());

	return 0;
}