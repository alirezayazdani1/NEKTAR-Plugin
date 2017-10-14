#include "vtkInformation.h"
#include "vtkMPIController.h"
#include "vtkNektarReader.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <iostream>

int main(int argc, char* argv[])
{
  if(argc != 2)
    {
    std::cerr << "Need to specify the nektar file location\n";
    return 1;
    }

  vtkMPIController* controller = vtkMPIController::New();
  controller->Initialize(&argc, &argv, 0);
  controller->SetGlobalController(controller);

  vtkNektarReader* reader = vtkNektarReader::New();
  //reader->SetFileName("/home/acbauer/DATA/DATA2KITWARE/test_cyl.nektar");
  reader->SetFileName(argv[1]);
  reader->SetUseProjection(1);
  reader->UpdateInformation();
#ifdef VTK5
  // this doesn't seem to actually apply the requested time step to the reader
  vtkInformation* outInfo = reader->GetOutputPortInformation(0);
  double time = 0.5;
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS(), &time, 1);
#else
  vtkInformation* outInfo = reader->GetOutputInformation(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), .5);
#endif
  reader->Update();
  reader->Delete();

  controller->Finalize(0);
  controller->SetGlobalController(0);
  controller->Delete();

  std::cerr << "finished\n";

  return 0;
}
