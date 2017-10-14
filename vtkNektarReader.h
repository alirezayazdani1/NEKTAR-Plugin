#ifndef __vtkNektarReader_h
#define __vtkNektarReader_h

#include <iostream>
#include <fstream>
#include <vector>

#include <nektar.h>
#include <gen_utils.h>

using namespace nektarTri;
using namespace nektarQuad;

#include "vtkUnstructuredGridAlgorithm.h"
class vtkPoints;
class vtkDataArraySelection;

#include "nektarObject.h"

class VTK_EXPORT vtkNektarReader : public vtkUnstructuredGridAlgorithm
{
 public:
  static vtkNektarReader *New();
  vtkTypeMacro(vtkNektarReader,vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  unsigned long GetMTime();

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  vtkSetStringMacro(DataFileName);
  vtkGetStringMacro(DataFileName);

  vtkGetMacro(NumberOfTimeSteps, int);
  // Description:
  // Returns the available range of valid integer time steps.
  vtkGetVector2Macro(TimeStepRange,int);
  vtkSetVector2Macro(TimeStepRange,int);

  // set/get whether projection should be used
  vtkSetMacro(UseProjection, int);
  vtkGetMacro(UseProjection, int);

 // set/get whether surface boundary should be extracted
  vtkSetMacro(ExtractBoundary, int);
  vtkGetMacro(ExtractBoundary, int);

   // set/get whether to use dynamic mesh (if is data is available)
  vtkSetMacro(DynamicMesh, int);
  vtkGetMacro(DynamicMesh, int);

  // set/get the scale to use for the deformation of the dynamic mesh
  void SetDynamicMeshScale(double);
  vtkGetMacro(DynamicMeshScale, double);
  
  // set/get the resolution to use for each element of the input grid
  void SetElementResolution(int);
  vtkGetMacro(ElementResolution, int);

  // set/get the resolution to use for each boundary element of the input grid (used for wss)
  void SetBoundaryResolution(int);
  vtkGetMacro(BoundaryResolution, int);

  // Description:
  // Get the number of point arrays available in the input.
  int GetNumberOfPointArrays(void);

  // Description:
  // Get the name of the  point array with the given index in
  // the input.
  const char* GetPointArrayName(int index);

  // Description:
  // Get/Set whether the point array with the given name or index is to
  // be read.
  int GetPointArrayStatus(const char* name);
  int GetPointArrayStatus(int index);
  void SetPointArrayStatus(const char* name, int status);

  // Description:
  // Turn on/off all point arrays.
  void DisableAllPointArrays();
  void EnableAllPointArrays();

  // Description:
  // Get the number of derived variables available in the input.
  int GetNumberOfDerivedVariableArrays(void);

  // Description:
  // Get the name of the  derived variable array with the given index in
  // the input.
  const char* GetDerivedVariableArrayName(int index);

  // Description:
  // Get/Set whether the derived variable array with the given name is to
  // be read.
  int GetDerivedVariableArrayStatus(const char* name);
  void SetDerivedVariableArrayStatus(const char* name, int status);

  // Description:
  // Turn on/off all derived variable arrays.
  void DisableAllDerivedVariableArrays();
  void EnableAllDerivedVariableArrays();

  // Description:
  // Get the names of any additional variables stored in the data
  int GetVariableNamesFromData();

 protected:
  vtkNektarReader();
  ~vtkNektarReader();

  char* FileName;
  char* DataFileName;
  int ElementResolution;
  int BoundaryResolution;
  int nfields;
  int my_patch_id;

  int num_vars;
  char** var_names;
  int num_der_vars;
  char** der_var_names;
  
  int* var_length;
  int* der_var_length;

  int dynamic_coord_index;
  int pressure_index;
  int velocity_index;
  int sm_displacement_index;
  int sm_acceleration_index;
  int sm_velocity_index;
  int num_used_scalars;
  int num_used_vectors;
  bool* use_field;

  //Tri* T;
  nektarList *myList;
  nektarObject *curObj;
  int displayed_step;
  int memory_step;
  int requested_step;

  FileList       fl;
  Element_List **master;
  Bndry         *Ubc;
  double        **WSS_all_vals;
  int            boundary_mem_step;

  static int next_patch_id;
  static bool NEED_TO_MANAGER_INIT;

  // to be read from .nektar file
  char paramFile[256];
  char reaFile[256];
  char rstFile[256];

  char  p_rea_file[256];
  char  p_rst_dir[256];
  char  p_rst_base[256];
  char  p_rst_ext[256];
  char  p_rst_file[256];
  char  p_rst_format[256];
  int   p_rst_start;
  int   p_rst_inc;
  int   p_rst_num;
  int   p_rst_digits;
  float p_time_start;
  float p_time_inc;

  void setActive();  // set my_patch_id as the active one
  static int getNextPatchID(){return(next_patch_id++);}

  vtkDataArraySelection* PointDataArraySelection;
  vtkDataArraySelection* DerivedVariableDataArraySelection;

  // update which fields from the data should be used, based on GUI
  void updateVariableStatus();
  // copy the data from nektar to pv
  void updateVtuData(vtkUnstructuredGrid* pv_ugrid, vtkUnstructuredGrid* pv_boundary_ugrid);
  void addCellsToContinuumMesh(int qa, double ***num);
  void interpolateAndCopyContinuumPoints(int alloc_res, int interp_res, vtkPoints* points);
  void interpolateAndCopyContinuumData(vtkUnstructuredGrid* pv_ugrid, double **data_array, int interp_res, int num_verts);
  void interpolateAndCopyBoundaryPoints(int alloc_res, int interp_res, vtkPoints* boundary_points);
  void interpolateAndCopyBoundaryData(int alloc_res, int num_verts, int interp_res);
  void addCellsToBoundaryMesh(int * boundary_index, int qa);
  void generateBoundaryConnectivity(int * boundary_index, int res);
  // see if the current object is missing data that was requested
  bool isObjectMissingData();
  // see if the current object matches the request
  bool objectMatchesRequest();
  // see if the current object has extra data than was requested
  bool objectHasExtraData();

  vtkUnstructuredGrid* UGrid;
  vtkUnstructuredGrid* Boundary_UGrid;
  bool CALC_GEOM_FLAG; // true = need to calculate continuum geometry; false = geom is up to date
  bool CALC_BOUNDARY_GEOM_FLAG; // true = need to calculate boundary geometry; false = boundary geom is up to date
  bool HAVE_BOUNDARY_GEOM_FLAG; // true = we have boundary geometry; false = geom has not been read yet

  bool READ_GEOM_FLAG; // true = need continuum geom from disk
  bool READ_BOUNDARY_GEOM_FLAG; // true = need boundary geom from disk
  bool IAM_INITIALLIZED;
  bool I_HAVE_DATA;
  bool FIRST_DATA;
  bool USE_MESH_ONLY;
  bool NEED_TO_BACKUP_MESH;

  //int TimeStep;
  int ActualTimeStep;
  int NumberOfTimeSteps;
  double TimeValue;
  int TimeStepRange[2];
  std::vector<double> TimeSteps;
  int UseProjection;
  int ExtractBoundary;
  int DynamicMesh;
  double DynamicMeshScale;

  // Time query function. Called by ExecuteInformation().
  // Fills the TimestepValues array.
  void GetAllTimes(vtkInformationVector*);

  // Description:
  // Populates the TIME_STEPS and TIME_RANGE keys based on file metadata.
  void AdvertiseTimeSteps( vtkInformation* outputInfo );

  virtual int RequestInformation(vtkInformation* request,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector);

  // Description:
  // This is called by the superclass.
  // This is the method you should override.
  virtual int RequestData(vtkInformation* request,
                          vtkInformationVector** inputVector,
                          vtkInformationVector* outputVector);

 private:
  vtkNektarReader(const vtkNektarReader&);  // Not implemented.
  void operator=(const vtkNektarReader&);  // Not implemented.
};

#endif
