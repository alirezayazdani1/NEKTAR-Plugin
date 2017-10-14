#include "vtkNektarReader.h"

#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkUnstructuredGrid.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTrivialProducer.h"
#include "vtkCleanUnstructuredGrid.h"
#include "vtkDataArraySelection.h"

#include "vtkFloatArray.h"
#include "vtkCellType.h"
#include "vtkPointData.h"
#include "vtkMultiProcessController.h"
#include "vtkNew.h"
#include "vtkSmartPointer.h"
#include "vtkTimerLog.h"

vtkStandardNewMacro(vtkNektarReader);

int  setup (FileList *f, Element_List **U, int *nftot, int Nsnapshots, bool mesh_only);
void ReadCopyField (FileList *f, Element_List **U);
void ReadAppendField (FileList *f, Element_List **U,  int start_field_index);
void Calc_StressTensor(Element_List **U, int displacement_index, int storage_index);
void Calc_Vort (FileList *f, Element_List **U, int nfields, int uvw_index, int p_index, int Snapshot_index);
void Calc_WSS(FileList *f, Element_List **E, Bndry *Ubc, int Snapshot_index, double** wss_vals=NULL);
int get_number_of_boundary_vertices(FileList *f, Element_List **E, Bndry *Ubc, int Snapshot_index);
void Set_Mesh(Element_List **U, int nfields, int iptr, double scale_val, int mode);
void Store_static_mesh(Element_List *U);

//----------------------------------------------------------------------------

int vtkNektarReader::next_patch_id = 0;

bool vtkNektarReader::NEED_TO_MANAGER_INIT = true;

vtkNektarReader::vtkNektarReader(){
  //this->DebugOn();
  vtkDebugMacro(<<"vtkNektarReader::vtkNektarReader(): ENTER");

  // by default assume filters have one input and one output
  // subclasses that deviate should modify this setting
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(2);

  this->FileName = 0;
  this->DataFileName = 0;
  this->ElementResolution = 1;
  this->BoundaryResolution = 1;
  this->my_patch_id = vtkNektarReader::getNextPatchID();
  vtkDebugMacro(<< "vtkNektarReader::vtkNektarReader(): my_patch_id = " << this->my_patch_id);

  this->UGrid = NULL;
  this->Boundary_UGrid = NULL;

  this->READ_GEOM_FLAG = true;
  this->READ_BOUNDARY_GEOM_FLAG = false;
  this->HAVE_BOUNDARY_GEOM_FLAG = false;
  this->CALC_GEOM_FLAG = true;
  this->CALC_BOUNDARY_GEOM_FLAG = false;
  this->IAM_INITIALLIZED = false;
  this->I_HAVE_DATA = false;
  this->FIRST_DATA = true;
  this->USE_MESH_ONLY = false;
  this->NEED_TO_BACKUP_MESH = false;

  this->UseProjection=0;
  this->ActualTimeStep = 0;
  this->DynamicMesh = 0;
  this->DynamicMeshScale = 1.0;
  this->TimeStepRange[0] = 0;
  this->TimeStepRange[1] = 0;
  this->NumberOfTimeSteps = 0;
  this->p_time_start = 0.0;
  this->p_time_inc = 0.0;
  this->displayed_step = -1;
  this->memory_step = -1;
  this->requested_step = -1;

  this->num_vars = 0;
  this->var_names = NULL;
  this->var_length = NULL;

  this->velocity_index = -1;
  this->pressure_index = -1;
  this->sm_displacement_index = -1;
  this->sm_acceleration_index = -1;
  this->sm_velocity_index = -1;
  
  this->dynamic_coord_index = -1;
  this->num_der_vars = 0;
  this->der_var_names = NULL;
  this->der_var_length = NULL;
  this->use_field = NULL;

 //this->nfields = 4;
  this->nfields = 0;
  memset (&this->fl, '\0', sizeof(FileList));

  this->fl.in.fp    =  stdin;
  this->fl.out.fp   =  stdout;
  this->fl.mesh.fp  =  stdin;

  this->fl.rea.name = new char[BUFSIZ];
  this->fl.in.name = new char[BUFSIZ];

  this->PointDataArraySelection = vtkDataArraySelection::New();
  //this->PointDataArraySelection->AddArray("Velocity");
  //this->PointDataArraySelection->AddArray("Pressure");

  this->DerivedVariableDataArraySelection = vtkDataArraySelection::New();
  //this->DerivedVariableDataArraySelection->AddArray("Vorticity");
  //this->DerivedVariableDataArraySelection->AddArray("lambda_2");
  //this->DerivedVariableDataArraySelection->AddArray("Wall Shear Stress");
  this->DisableAllDerivedVariableArrays();

  //this->master = (Element_List **) malloc((2)*this->nfields*sizeof(Element_List *));
  this->master = NULL;
  this->Ubc = NULL;
  this->WSS_all_vals = NULL;
  this->boundary_mem_step = -1;

  this->myList = new nektarList();

  vtkDebugMacro(<<"vtkNektarReader::vtkNektarReader(): EXIT");
}

//----------------------------------------------------------------------------
vtkNektarReader::~vtkNektarReader()
{
  vtkDebugMacro(<<"vtkNektarReader::~vtkNektarReader(): ENTER");

  int i;
  this->SetFileName(0);
  this->SetDataFileName(0);
  this->PointDataArraySelection->Delete();
  this->DerivedVariableDataArraySelection->Delete();

  if(this->UGrid)
    {
    vtkDebugMacro(<<"vtkNektarReader::~vtkNektarReader(): ugrid not null, delete it");
    this->UGrid->Delete();
    }
  if(this->Boundary_UGrid)
    {
    this->Boundary_UGrid->Delete();
    }
  if(this->WSS_all_vals)
    {
    delete this->WSS_all_vals;
    }

  if (this->myList)
    {
    myList->~nektarList();
    }

  if(this->num_vars>0)
  {
      free(this->var_length);
      for (i=0; i<this->num_vars; i++)
      {
	  free(this->var_names[i]);
      }
      free(this->var_names);
  }
  if(this->num_der_vars>0)
  {
      free(this->der_var_length);
      for (i=0; i<this->num_der_vars; i++)
      {
	  free(this->der_var_names[i]);
      }
      free(this->der_var_names);
      this->der_var_names = NULL;
  }
  if(this->use_field)
  {
      free(this->use_field);
      this->use_field= NULL;
  }

  vtkDebugMacro(<<"vtkNektarReader::~vtkNektarReader(): clean up master");
  /// We need to clean up master, Ubc, WSS_all_vals, others?
	  
// jai -- this is not working, need to investigate
	  
//   if(this->master)
//   {
// 	  for(i=0; i<this->nfields; i++)
// 	  {
// 		  this->master[i]->Element_List::~Element_List();
// 		  vtkDebugMacro(<<"vtkNektarReader::~vtkNektarReader(): master ["<<i<<"]->~Element_List() called");
// 	  }
// 	  free(this->master);
// 	  vtkDebugMacro(<<"vtkNektarReader::~vtkNektarReader(): master freed");
//   }

  vtkDebugMacro(<<"vtkNektarReader::~vtkNektarReader(): EXIT");
}

//----------------------------------------------------------------------------
void vtkNektarReader::setActive()
{
  iparam_set("IDpatch", this->my_patch_id);
}

//----------------------------------------------------------------------------
void vtkNektarReader::GetAllTimes(vtkInformationVector *outputVector)
{
  FILE* dfPtr = NULL;
  char dfName[265];
  char* scan_ret;
  char* p;
  char* p2;
  char param[32];
  char paramLine[256];
  int file_index;
  float test_time_val;

  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkInformation* outInfo1 = outputVector->GetInformationObject(1);

  this->TimeStepRange[0] = 0;
  this->TimeStepRange[1] = this->NumberOfTimeSteps-1;

  vtkDebugMacro(<< "vtkNektarReader::GetAllTimes: this->NumberOfTimeSteps = "<<
                this->NumberOfTimeSteps);

  this->TimeSteps.resize(this->NumberOfTimeSteps);

  if(!this->USE_MESH_ONLY)
    {
    if(this->p_time_inc == 0.0)
      {
      for (int i=0; i<(this->NumberOfTimeSteps); i++)
        {
        file_index = this->p_rst_start + (this->p_rst_inc*i);
        sprintf(dfName, this->p_rst_format, file_index );
        dfPtr = fopen(dfName, "r");
        if(!dfPtr)
          {
          vtkErrorMacro(<< "Failed to open file: "<< dfName);
          }
        // skip the first 5 lines
        for(int j=0; j<5; j++)
          {
          scan_ret = fgets(paramLine, 256, dfPtr);
          }
        scan_ret = fgets(paramLine, 256, dfPtr);

        //fprintf(stderr, "Line: \'%s\'\n", paramLine);
        //sscanf(paramLine, "%f", &(this->TimeSteps[i]));
        sscanf(paramLine, "%f", &test_time_val);
        this->TimeSteps[i] = test_time_val;

        fclose(dfPtr);
        dfPtr = NULL;
        //fprintf(stderr, "File: %s : time: %f\n", dfName, this->TimeSteps[i]);
        }
      }
    else
      {
      double cur_time = this->p_time_start;
      for (int i=0; i<(this->NumberOfTimeSteps); i++)
        {
        this->TimeSteps[i] = cur_time;
        cur_time += this->p_time_inc;
        }

      }
    }  // from if(!this->USE_MESH_ONLY)
  else
    {
    this->TimeSteps[0] = 0.0;
    }
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
               &(this->TimeSteps[0]),
               this->NumberOfTimeSteps);

  outInfo1->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
                &(this->TimeSteps[0]),
                this->NumberOfTimeSteps);

  double timeRange[2];
  timeRange[0] = this->TimeSteps[0];
  timeRange[1] = this->TimeSteps[this->NumberOfTimeSteps-1];

  vtkDebugMacro(<< "vtkNektarReader::GetAllTimes: timeRange[0] = "<<timeRange[0]<< ", timeRange[1] = "<< timeRange[1]);

  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),
               timeRange, 2);
  outInfo1->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),
                timeRange, 2);
} // vtkNektarReader::GetAllTimes()

//----------------------------------------------------------------------------
unsigned long vtkNektarReader::GetMTime()
{
  unsigned long mTime = this->Superclass::GetMTime();
  unsigned long time;

  time = this->PointDataArraySelection->GetMTime();
  mTime = ( time > mTime ? time : mTime );

  time = this->DerivedVariableDataArraySelection->GetMTime();
  mTime = ( time > mTime ? time : mTime );

  return mTime;
}

//----------------------------------------------------------------------------
void vtkNektarReader::SetDynamicMeshScale(double val)
{
     
    this->DynamicMeshScale = val;
    //iparam_set("VIS_RES", ElementResolution);

    // *** is this still always true?  I'm thinking not, but it
    //     may not matter....
    this->CALC_GEOM_FLAG = true;
    this->CALC_BOUNDARY_GEOM_FLAG = true;
    this->Modified();
    vtkDebugMacro(<<"vtkNektarReader::SetMeshScale: DynamicMeshScale= "<< this->DynamicMeshScale);
    
}


//----------------------------------------------------------------------------
void vtkNektarReader::SetElementResolution(int val)
{
  if(val >0 && val<11)
    {
    this->ElementResolution = val;
    //iparam_set("VIS_RES", ElementResolution);

    // *** is this still always true?  I'm thinking not, but it
    //     may not matter....
    this->CALC_GEOM_FLAG = true;
    this->Modified();
    vtkDebugMacro(<<"vtkNektarReader::SetElementResolution: CALC_GEOM_FLAG now true (need to calculate) , ElementResolution= "<< this->ElementResolution);
    }
}

//----------------------------------------------------------------------------
void vtkNektarReader::SetBoundaryResolution(int val)
{
  if(val >0 && val<11)
    {
    this->BoundaryResolution = val;
    //iparam_set("VIS_RES", ElementResolution);

    // *** is this still always true?  I'm thinking not, but it
    //     may not matter....
    this->CALC_BOUNDARY_GEOM_FLAG = true;
    this->Modified();
    vtkDebugMacro(<<"vtkNektarReader::SetBoundaryResolution: CALC_BOUNDARY_GEOM_FLAG now true (need to calculate) , BoundaryResolution= "<< this->BoundaryResolution);
    }
}

//----------------------------------------------------------------------------
void vtkNektarReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
int vtkNektarReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}

//----------------------------------------------------------------------------
const char* vtkNektarReader::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}

//----------------------------------------------------------------------------
int vtkNektarReader::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}

//----------------------------------------------------------------------------
int vtkNektarReader::GetPointArrayStatus(int index)
{
  return this->PointDataArraySelection->GetArraySetting(index);
}

//----------------------------------------------------------------------------
void vtkNektarReader::SetPointArrayStatus(const char* name, int status)
{
  if(status)
    {
    this->PointDataArraySelection->EnableArray(name);
    }
  else
    {
    this->PointDataArraySelection->DisableArray(name);
    }
}

//----------------------------------------------------------------------------
void vtkNektarReader::EnableAllPointArrays()
{
  this->PointDataArraySelection->EnableAllArrays();
}

//----------------------------------------------------------------------------
void vtkNektarReader::DisableAllPointArrays()
{
  this->PointDataArraySelection->DisableAllArrays();
}

//----------------------------------------------------------------------------
int vtkNektarReader::GetNumberOfDerivedVariableArrays()
{
  return this->DerivedVariableDataArraySelection->GetNumberOfArrays();
}

//----------------------------------------------------------------------------
const char* vtkNektarReader::GetDerivedVariableArrayName(int index)
{
  return this->DerivedVariableDataArraySelection->GetArrayName(index);
}

//----------------------------------------------------------------------------
int vtkNektarReader::GetDerivedVariableArrayStatus(const char* name)
{
  return this->DerivedVariableDataArraySelection->ArrayIsEnabled(name);
}

//----------------------------------------------------------------------------
void vtkNektarReader::SetDerivedVariableArrayStatus(const char* name, int status)
{
  if(status)
    {
    this->DerivedVariableDataArraySelection->EnableArray(name);
    }
  else
    {
    this->DerivedVariableDataArraySelection->DisableArray(name);
    }
}

//----------------------------------------------------------------------------
void vtkNektarReader::EnableAllDerivedVariableArrays()
{
  this->DerivedVariableDataArraySelection->EnableAllArrays();
}

//----------------------------------------------------------------------------
void vtkNektarReader::DisableAllDerivedVariableArrays()
{
  this->DerivedVariableDataArraySelection->DisableAllArrays();
}

//----------------------------------------------------------------------------

void vtkNektarReader::updateVariableStatus()
{
    int id = 0;
    int i;
    int my_rank = vtkMultiProcessController::GetGlobalController()->GetLocalProcessId();
    //vtkDebugMacro(<<"vtkNektarReader::updateVariableStatus: Rank: "<<my_rank<< " ENTER ");
    this->num_used_vectors=0;
    this->num_used_scalars=0;

    // set all fields to false
    for(i=0; i<this->nfields; i++)
    {
	this->use_field[i] = false;
    }

    // if a variable is used, set it to true
    for(i=0; i<this->num_vars; i++)
    {
	if(this->GetPointArrayStatus(i) == 1)
	{
	    this->use_field[id] = true;
	    // if it is a vector, set the next two components to true as well
	    // also increment the number of vectors or scalars used, accordingly
	    if(this->var_length[i] == 3)
	    {
		this->use_field[id+1] = true;
		this->use_field[id+2] = true;
		this->num_used_vectors++;
	    }
	    else
	    {
		this->num_used_scalars++;
	    }
	}
	if(this->var_length[i] == 1)
	{
	    id++;
	}
	else if(this->var_length[i] == 3)
	{
	    id+=3;
	}
    }
//     for(i=0; i<this->nfields; i++)
//     {
// 	vtkDebugMacro(<<"vtkNektarReader::updateVariableStatus: Rank: "<<my_rank<< ": this->use_field["<<i<<"]: "<<this->use_field[i]);
//     }
    vtkDebugMacro(<<"vtkNektarReader::updateVariableStatus: Rank: "<<my_rank<< ": this->num_used_scalars= "<<this->num_used_scalars<<" : this->num_used_vectors= "<<this->num_used_vectors);
}

//----------------------------------------------------------------------------
int vtkNektarReader::GetVariableNamesFromData()
{
  FILE* dfPtr = NULL;
  char  dfName[265];
  char* scan_ret;
  char  paramLine[256];
  int   ind = 0;
  char  l_var_name[2];
  char  l_duplicate_var_name[8];
  int   dup_var_number=1;

  bool this_var_added = false;
  
  this->num_vars = 0;

  l_var_name[1] = '\0';
  int my_rank = vtkMultiProcessController::GetGlobalController()->GetLocalProcessId();
  if(my_rank == 0)
    {
    sprintf(dfName, this->p_rst_format, this->p_rst_start );
    dfPtr = fopen(dfName, "r");
    if(!dfPtr) vtkErrorMacro(<< "Failed to open file: "<< dfName);
    // skip the first 8 lines
    for(int j=0; j<8; j++)
      {
      scan_ret = fgets(paramLine, 256, dfPtr);
      }
    scan_ret = fgets(paramLine, 256, dfPtr);

    fclose(dfPtr);
    dfPtr = NULL;

    vtkDebugMacro(<< "vtkNektarReader::GetVariableNamesFromData:before bcast my_rank: "<<my_rank<< "  paramLine = "<< paramLine);
    }

  vtkMultiProcessController::GetGlobalController()->Broadcast(paramLine, 256, 0);
  vtkDebugMacro(<< "vtkNektarReader::GetVariableNamesFromData:after  bcast my_rank: "<<my_rank<< "  paramLine = "<< paramLine);

  char* space_index = index(paramLine, ' ');
  *space_index = '\0';
  int len = strlen(paramLine);
  vtkDebugMacro(<< "vtkNektarReader::GetVariableNamesFromData:after strlen my_rank: "<<my_rank<< "  paramLine = \'"<< paramLine<< "\'  len= "<< len);

  // allocate space for variable names and lengths, 
  // will be at most one for each character in paramLine (if all are scalars)
  this->var_names =  (char**) malloc(len* sizeof(char*));
  this->var_length =  (int*) malloc(len* sizeof(int));

  // if the variable list contains a 'p' (pressure), then this is a fluid data set
  if(strchr(paramLine, 'p')!= NULL)
  {
      
      while(ind<len)
      {
          this_var_added = false;
	  vtkDebugMacro(<< "vtkNektarReader::GetVariableNamesFromData: my_rank: "<<my_rank<< "  paramLine["<<ind<<"] = "<< paramLine[ind]<<" : this->num_vars = "<< this->num_vars);
	  switch(paramLine[ind])
	  {

	      case 'p':
		  this->PointDataArraySelection->AddArray("Pressure");
		  this->var_names[this->num_vars] = strdup("Pressure");
		  this->var_length[this->num_vars] = 1; // this is a scalar
		  this->pressure_index=ind;
		  ind++;
		  this_var_added = true;
		  break;
	      case 'u':
		  if((ind+2 < len) && (paramLine[ind+1] == 'v' && paramLine[ind+2] == 'w'))
		  {
		      this->PointDataArraySelection->AddArray("Fluid Velocity");
		      this->var_names[this->num_vars] = strdup("Fluid Velocity");
		      this->var_length[this->num_vars] = 3; // this is a vector
		      this->velocity_index = ind;
		      ind+=3;
		      this_var_added = true;

		      // if we have Velocity, we can derive these quantities
		      this->num_der_vars = 3;
		      this->der_var_names =  (char**) malloc(this->num_der_vars * sizeof(char*));
		      this->der_var_length =  (int*) malloc(this->num_der_vars * sizeof(int));

		      this->DerivedVariableDataArraySelection->AddArray("Vorticity");
		      this->der_var_names[0] = strdup("Vorticity");
		      this->der_var_length[0] = 3; // this is a vector
		      this->DerivedVariableDataArraySelection->AddArray("lambda_2");
		      this->der_var_names[1] = strdup("lambda_2");
		      this->der_var_length[1] = 1; // this is a scalar
		      this->DerivedVariableDataArraySelection->AddArray("Wall Shear Stress");
		      this->der_var_names[2] = strdup("Wall Shear Stress");
		      this->der_var_length[2] = 1; // this is a scalar
		      
		      //this->DisableAllDerivedVariableArrays();
		  }
		  break;
	      case 'i':
		  if((ind+2 < len) && (paramLine[ind+1] == 'j' && paramLine[ind+2] == 'k'))
		  {
		      this->PointDataArraySelection->AddArray("Fluid Mesh Displacement");
		      this->var_names[this->num_vars] = strdup("Fluid Mesh Displacement");
		      this->var_length[this->num_vars] = 3; // this is a vector
		      this->dynamic_coord_index = ind;
		      this->NEED_TO_BACKUP_MESH = true;
		      ind+=3;
		      this_var_added = true;
		  }
		  break;
	      case 'r':
		  if((ind+2 < len) && (paramLine[ind+1] == 's' && paramLine[ind+2] == 't'))
		  {
		      this->PointDataArraySelection->AddArray("Fluid Mesh Velocity");
		      this->var_names[this->num_vars] = strdup("Fluid Mesh Velocity");
		      this->var_length[this->num_vars] = 3; // this is a vector
		      ind+=3;
		      this_var_added = true;
		  }
		  break;
	      default:
		  l_var_name[0] = paramLine[ind];
// 		  this->PointDataArraySelection->AddArray(l_var_name);
// 		  this->var_names[this->num_vars] = strdup(l_var_name);
// 		  this->var_length[this->num_vars] = 1; // this is a scalar
// 		  ind++;
		  break;
	  }
	  if (!this_var_added)
	  {
		  l_var_name[0] = paramLine[ind];
		  if (this->PointDataArraySelection->ArrayExists(l_var_name))
		  {
		    sprintf(l_duplicate_var_name, "%s_%02d", l_var_name, dup_var_number++);
		    this->PointDataArraySelection->AddArray(l_duplicate_var_name);
		    this->var_names[this->num_vars] = strdup(l_duplicate_var_name);
		  }
		  else
		  {
		    this->PointDataArraySelection->AddArray(l_var_name);
		    this->var_names[this->num_vars] = strdup(l_var_name);
		  }
		  this->var_length[this->num_vars] = 1; // this is a scalar
		  ind++;		  
	  }
	  this->num_vars++;
      }
      /*
      this->num_der_vars = 3;
      this->der_var_names =  (char**) malloc(this->num_der_vars * sizeof(char*));
      this->der_var_length =  (int*) malloc(this->num_der_vars * sizeof(int));


      this->DerivedVariableDataArraySelection->AddArray("Vorticity");
      this->der_var_names[0] = strdup("Vorticity");
      this->der_var_length[0] = 3; // this is a vector
      this->DerivedVariableDataArraySelection->AddArray("lambda_2");
      this->der_var_names[1] = strdup("lambda_2");
      this->der_var_length[1] = 1; // this is a scalar
      this->DerivedVariableDataArraySelection->AddArray("Wall Shear Stress");
      this->der_var_names[2] = strdup("Wall Shear Stress");
      this->der_var_length[2] = 1; // this is a scalar

      //this->DisableAllDerivedVariableArrays();
      */

  }
  else  /// this is a solid (or, at least not a fluid...)
  {
      int add_stress_tensor = 0;
      while(ind<len)
      {
          this_var_added = false;
	  vtkDebugMacro(<< "vtkNektarReader::GetVariableNamesFromData: my_rank: "<<my_rank<< "  paramLine["<<ind<<"] = "<< paramLine[ind]<<" : this->num_vars = "<< this->num_vars);
	  switch(paramLine[ind])
	  {
	      case 'u':
		  if((ind+2 < len) && (paramLine[ind+1] == 'v' && paramLine[ind+2] == 'w'))
		  {
		      this->PointDataArraySelection->AddArray("Solid Mesh Displacement");
		      this->var_names[this->num_vars] = strdup("Solid Mesh Displacement");
		      this->var_length[this->num_vars] = 3; // this is a vector
		      this->dynamic_coord_index = ind;
		      this->sm_displacement_index = ind;
		      this->NEED_TO_BACKUP_MESH = true;
		      ind+=3;
		      add_stress_tensor++;
		      this_var_added = true;
		  }
		  break;
	      case 'i':
		  if((ind+2 < len) && (paramLine[ind+1] == 'j' && paramLine[ind+2] == 'k'))
		  {
		      this->PointDataArraySelection->AddArray("Solid Mesh Acceleration");
		      this->var_names[this->num_vars] = strdup("Solid Mesh Acceleration");
		      this->var_length[this->num_vars] = 3; // this is a vector
		      this->sm_acceleration_index = ind;
		      ind+=3;
		      this_var_added = true;		      
		  }
		  break;
	      case 'r':
		  if((ind+2 < len) && (paramLine[ind+1] == 's' && paramLine[ind+2] == 't'))
		  {
		      this->PointDataArraySelection->AddArray("Solid Mesh Velocity");
		      this->var_names[this->num_vars] = strdup("Solid Mesh Velocity");
		      this->var_length[this->num_vars] = 3; // this is a vector

		      this->sm_velocity_index = ind;
		      add_stress_tensor++;
		      ind+=3;
		      this_var_added = true;
		  }
		  break;
	      default:
// - moved this down to catch possiblity that above variable 'u' may not be followed by 'v' and 'w',
//     same for 'i' and 'r'		      
// 		  l_var_name[0] = paramLine[ind];
// 		  this->PointDataArraySelection->AddArray(l_var_name);
// 		  this->var_names[this->num_vars] = strdup(l_var_name);
// 		  this->var_length[this->num_vars] = 1; // this is a scalar
// 		  ind++;
		  break;
	  }
	  if (!this_var_added)
	  {
		  l_var_name[0] = paramLine[ind];
		  this->PointDataArraySelection->AddArray(l_var_name);
		  this->var_names[this->num_vars] = strdup(l_var_name);
		  this->var_length[this->num_vars] = 1; // this is a scalar
		  ind++;		  
	  }
	  this->num_vars++;
      }
      if(2 == add_stress_tensor)
      {
	      this->num_der_vars = 1;
	      this->der_var_names =  (char**) malloc(this->num_der_vars * sizeof(char*));
	      this->der_var_length =  (int*) malloc(this->num_der_vars * sizeof(int));
	      
	      this->DerivedVariableDataArraySelection->AddArray("Stress Tensor");
	      this->der_var_names[0] = strdup("Stress Tensor");
	      this->der_var_length[0] = 3; // this is a vector

	      //this->DisableAllDerivedVariableArrays();
      }
  }

  return len;
}

//----------------------------------------------------------------------------
int vtkNektarReader::RequestInformation(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector),
  vtkInformationVector* outputVector)
{
  vtkDebugMacro(<<"vtkNektarReader::RequestInformation(): ENTER");

  int numArrays;
  int nprocs;
  int mytid;
  FILE* inPtr = NULL;
  char* scan_ret;
  char* p;
  char* p2;
  char param[32];
  char paramLine[256];
  double timer_diff;


  if(!this->IAM_INITIALLIZED)
    {
    //print the name of the file we're supposed to open
    vtkDebugMacro(<< "vtkNektarReader::RequestInformation: FileName: " << this->GetFileName());

    inPtr = fopen(this->GetFileName(), "r");
    if(inPtr == NULL)
      {
      vtkErrorMacro(<< "vtkNektarReader::RequestInformation: could not open file: "<< this->GetFileName());
      return(0);
      }

    this->p_rst_digits = 0;

    scan_ret = fgets(paramLine, 256, inPtr);

    while(scan_ret != NULL)
      {
      p=strchr(paramLine, ':');

      while(p != NULL && !this->USE_MESH_ONLY)
        {
        *p = '\0';
        sscanf(paramLine, "%s", param);
        p = p+2;

        p2=strchr(p, '\n');
        *p2= '\0';

        if(strcasecmp(param, "REA_FILE") == 0)
          {
          while (*p == ' ')
            {
            p++;
            }
          strcpy(this->p_rea_file, p);
          }
        else if(strcasecmp(param, "RST_DIR") == 0)
          {
          while (*p == ' ')
            {
            p++;
            }
          strcpy(this->p_rst_dir, p);
          if(strcasecmp(this->p_rst_dir, "NULL") == 0)
            {
            this->USE_MESH_ONLY = true;
            this->p_rst_num = 1;
	    this->SetExtractBoundary(1);
	    vtkDebugMacro(<< "vtkNektarReader::RequestInformation: this->GetExtractBoundary(): " << this->GetExtractBoundary());

	    if( !this->HAVE_BOUNDARY_GEOM_FLAG)
	    {
		this->READ_BOUNDARY_GEOM_FLAG = true;
	    }

            p= NULL;
            //fprintf(stderr, "found RST_DIR: %s\n", this->p_rst_dir);
            fflush(stderr);
            break;
            }
          }
        else if(strcasecmp(param, "RST_BASE") == 0)
          {
          while (*p == ' ')
            {
            p++;
            }
          strcpy(this->p_rst_base, p);
          }
        else if(strcasecmp(param, "RST_EXT") == 0)
          {
          while (*p == ' ')
            {
            p++;
            }
          strcpy(this->p_rst_ext, p);
          }
        else if(strcasecmp(param, "RST_START") == 0)
          {
          sscanf(p, "%d", &this->p_rst_start);
          }
        else if(strcasecmp(param, "RST_INC") == 0)
          {
          sscanf(p, "%d", &this->p_rst_inc);
          }
        else if(strcasecmp(param, "RST_NUM") == 0)
          {
          sscanf(p, "%d", &this->p_rst_num);
          }
        else if(strcasecmp(param, "RST_DIGITS") == 0)
          {
          sscanf(p, "%d", &this->p_rst_digits);
          }
        else if(strcasecmp(param, "TIME_START") == 0)
          {
          sscanf(p, "%f", &this->p_time_start);
          }
        else if(strcasecmp(param, "TIME_INC") == 0)
          {
          sscanf(p, "%f", &this->p_time_inc);
          }
        scan_ret = fgets(paramLine, 256, inPtr);
        if(scan_ret)
          {
          p=strchr(paramLine, ':');
          }
        else
          {
          p= NULL;
          }
        } // while(p!=NULL && !this->USE_MESH_ONLY)

      scan_ret = fgets(paramLine, 256, inPtr);
      }// while(ret != NULL)

    //this->p_rst_current = this->p_rst_start;

    vtkDebugMacro(<<"REA_FILE:   \'"<<this->p_rea_file<<"\'");
    vtkDebugMacro(<<"RST_DIR:    \'"<<this->p_rst_dir<<"\'");
    if(!this->USE_MESH_ONLY)
      {
      vtkDebugMacro(<<"RST_BASE:   \'"<<this->p_rst_base<<"\'");
      vtkDebugMacro(<<"RST_EXT:    \'"<<this->p_rst_ext<<"\'");
      vtkDebugMacro(<<"RST_START:  "<<this->p_rst_start);
      vtkDebugMacro(<<"RST_INC:    "<<this->p_rst_inc);
      vtkDebugMacro(<<"RST_NUM:    "<<this->p_rst_num);
      sprintf(this->p_rst_format, "%s/%s%%%dd.%s",
              this->p_rst_dir,
              this->p_rst_base,
              this->p_rst_digits,
              this->p_rst_ext);
      vtkDebugMacro(<<"RST_FORMAT: \'"<<this->p_rst_format<<"\'");
      }
    this->NumberOfTimeSteps = this->p_rst_num;

    fclose(inPtr);

    vtkNew<vtkTimerLog> timer;
    timer->StartTimer();
    this->GetAllTimes(outputVector);
    timer->StopTimer();
    timer_diff = timer->GetElapsedTime();

    if (!this->USE_MESH_ONLY)
      {
	  this->nfields = this->GetVariableNamesFromData();
      }

    this->master = (Element_List **) malloc(this->nfields*sizeof(Element_List *));
    this->use_field = (bool*) malloc(this->nfields*sizeof(bool));
    vtkDebugMacro(<<"Rank: "<<vtkMultiProcessController::GetGlobalController()->GetLocalProcessId()
                  <<" :: GetAllTimeSteps time: "<< timer_diff);

//      sprintf(this->fl.rea.name,"%s", this->reaFile);
//      this->SetDataFileName(this->reaFile);
    sprintf(this->fl.rea.name,"%s", this->p_rea_file);
    this->SetDataFileName(this->p_rea_file);

    this->fl.rea.fp = fopen(this->fl.rea.name,"r");
    if (!this->fl.rea.fp)
      {
      error_msg(Restart: no REA file read);
      }

    if(!this->USE_MESH_ONLY)
      {
      sprintf(this->p_rst_file, this->p_rst_format, this->p_rst_start);
      sprintf(this->fl.in.name,"%s", this->p_rst_file);
      this->fl.in.fp = fopen(this->fl.in.name,"r");
      if (!this->fl.in.fp)
        {
        error_msg(Restart: no dumps read from restart file);
        }
      }

    vtkInformation *outInfo0 =
      outputVector->GetInformationObject(0);
    //outInfo0->Set(vtkStreamingDemandDrivenPipeline::
    //              MAXIMUM_NUMBER_OF_PIECES(), -1);
    outInfo0->Set(vtkAlgorithm::CAN_HANDLE_PIECE_REQUEST(), 1);

    vtkInformation *outInfo1 =
      outputVector->GetInformationObject(1);
    //outInfo1->Set(vtkStreamingDemandDrivenPipeline::
    //              MAXIMUM_NUMBER_OF_PIECES(), -1);
    outInfo1->Set(vtkAlgorithm::CAN_HANDLE_PIECE_REQUEST(), 1);

    this->IAM_INITIALLIZED = true;
    }// if(!this->IAM_INITIALLIZED)

  vtkDebugMacro(<<"vtkNektarReader::RequestInformation(): EXIT");

  return 1;
}

//----------------------------------------------------------------------------
// This is the superclasses style of Execute method.  Convert it into
// an imaging style Execute method.
int vtkNektarReader::RequestData(
  vtkInformation* request,
  vtkInformationVector** vtkNotUsed( inputVector ),
  vtkInformationVector* outputVector)
{
  double timer_diff;
  double total_timer_diff;

  vtkNew<vtkTimerLog> total_timer;
  total_timer->StartTimer();
  // the default implimentation is to do what the old pipeline did find what
  // output is requesting the data, and pass that into ExecuteData

  // which output port did the request come from
  int outputPort =
    request->Get(vtkDemandDrivenPipeline::FROM_OUTPUT_PORT());

  vtkDebugMacro(<<"RequestData: ENTER: outputPort = "<< outputPort);
  
  // if output port is negative then that means this filter is calling the
  // update directly, in that case just assume port 0
  if (outputPort == -1)
    {
    outputPort = 0;
    }

  // get the data object
  vtkInformation *outInfo =
    outputVector->GetInformationObject(0);     //(outputPort);

  vtkInformation *outInfo1 =
    outputVector->GetInformationObject(1);

  vtkInformation *outInfoArray[2];
  outInfoArray[0] = outInfo;
  outInfoArray[1] = outInfo1;
	  
  vtkInformation *requesterInfo =
    outputVector->GetInformationObject(outputPort);

  int tsLength =
    requesterInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

  double* steps =
    requesterInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

  vtkDebugMacro(<<"RequestData: tsLength= "<< tsLength);

  // Update the status of the requested variables
  this->updateVariableStatus();

  double l_time_val_0 = 0.0;
  double l_time_val_1 = 0.0;
  
  // Check if a particular time was requested.
  bool hasTimeValue = false;
#ifdef VTK5
  vtkDebugMacro(<<"RequestData: #ifdef VTK5 ");
  if(requesterInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS()))
    {
    // Get the requested time step. We only supprt requests of a single time
    // step in this reader right now
    double *requestedTimeSteps =
      requesterInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEPS());
    this->TimeValue = requestedTimeSteps[0];
    hasTimeValue = true;
    }
#else
  vtkDebugMacro(<<"RequestData: #else (!ifdef VTK5) ");
  // Collect the time step requested
  vtkInformationDoubleKey* timeKey =
    vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP();

  if (outInfoArray[outputPort]->Has(timeKey))
    {
    this->TimeValue = outInfoArray[outputPort]->Get(timeKey);
    hasTimeValue = true;
    }
#endif
  if(hasTimeValue)
    {
    vtkDebugMacro(<<"RequestData: this->TimeValue= "<< this->TimeValue);

    //find the timestep with the closest value to the requested time value
    int closestStep=0;
    double minDist=-1;
    for (int cnt=0;cnt<tsLength;cnt++)
      {
      //fprintf(stderr, "RequestData: steps[%d]=%f\n", cnt, steps[cnt]);
      double tdist=(steps[cnt]-this->TimeValue>this->TimeValue-steps[cnt])?steps[cnt]-this->TimeValue:this->TimeValue-steps[cnt];
      if (minDist<0 || tdist<minDist)
        {
        minDist=tdist;
        closestStep=cnt;
        }
      }
    this->ActualTimeStep=closestStep;
    }

  vtkDebugMacro(<<"RequestData: this->ActualTimeStep= "<< this->ActualTimeStep);

  // Force TimeStep into the "known good" range. Although this
  if ( this->ActualTimeStep < this->TimeStepRange[0] )
    {
    this->ActualTimeStep = this->TimeStepRange[0];
    }
  else if ( this->ActualTimeStep > this->TimeStepRange[1] )
    {
    this->ActualTimeStep = this->TimeStepRange[1];
    }

  int my_rank = vtkMultiProcessController::GetGlobalController()->GetLocalProcessId();
  vtkDebugMacro(<<"RequestData: ENTER: rank: "<< my_rank << "  outputPort: "
                << outputPort << "  this->ActualTimeStep = "<< this->ActualTimeStep);

  // if the user has requested wss, and we have not read the boundary geometry
  // before, set READ_BOUNDARY_GEOM_FLAG to true.  This will only happen once.
  if(this->GetDerivedVariableArrayStatus("Wall Shear Stress") &&
     !this->HAVE_BOUNDARY_GEOM_FLAG)
    {
    this->READ_BOUNDARY_GEOM_FLAG = true;
    }

  // if user said to extract boundary, and we haven't yet, set READ_BOUNDARY_GEOM_FLAG to true.  This will only happen once.
  if( this->GetExtractBoundary() && !this->HAVE_BOUNDARY_GEOM_FLAG)
  {
      this->READ_BOUNDARY_GEOM_FLAG = true;
  }

#if 0
  /* Paris: Test what happens if we read boundary geometry anyways */
  if( !this->HAVE_BOUNDARY_GEOM_FLAG)
  {
      this->READ_BOUNDARY_GEOM_FLAG = true;
  }
#endif

  vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid* boundary_ugrid = vtkUnstructuredGrid::SafeDownCast(outInfo1->Get(vtkDataObject::DATA_OBJECT()));

  // Save the time value in the output (ugrid) data information.
  if (steps)
    {
#ifdef VTK5
    ugrid->GetInformation()->Set(vtkDataObject::DATA_TIME_STEPS(),
                                 steps+this->ActualTimeStep, 1);
    boundary_ugrid->GetInformation()->Set(vtkDataObject::DATA_TIME_STEPS(),
                                     steps+this->ActualTimeStep, 1);
#else
    ugrid->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(),
                                 steps[this->ActualTimeStep]);
    boundary_ugrid->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(),
                                     steps[this->ActualTimeStep]);
#endif
    }

//  int new_rst_val = this->p_rst_start + (this->p_rst_inc* this->ActualTimeStep);
  this->requested_step = this->p_rst_start + (this->p_rst_inc* this->ActualTimeStep);

  //  if the step being displayed is different than the one requested
  //if(this->displayed_step != this->requested_step)
    {
    // get the requested object from the list, if the ugrid in the object is NULL
    // then we have not loaded it yet
    this->curObj = this->myList->getObject(this->requested_step);

    // need to compare the request with the curObj
    // if the resolutions are not the same, or any variable in the request is true,
    // but false in the object, then we will need the data to be in memory
    
 //    if( (this->curObj->resolution != this->ElementResolution) ||
//         (this->curObj->boundary_resolution != this->BoundaryResolution) ||
//         (this->curObj->use_projection != this->UseProjection) ||
//         (this->GetPointArrayStatus("Pressure") && !this->curObj->pressure) ||
//         (this->GetPointArrayStatus("Velocity") && !this->curObj->velocity) ||
//         (this->GetDerivedVariableArrayStatus("Vorticity") && !this->curObj->vorticity) ||
//         (this->GetDerivedVariableArrayStatus("lambda_2") && !this->curObj->lambda_2) ||
//         (this->GetDerivedVariableArrayStatus("Wall Shear Stress") && !this->curObj->wss) )

    if(this->isObjectMissingData())
      {
      // if the step in memory is different than the step requested
      if(this->requested_step != this->memory_step)
        {
        this->I_HAVE_DATA = false;
        }
      // if the current object uses dynamic mesh and the request does not (or vice versa), will need to recalculate the geometry
      vtkDebugMacro(<< "vtkNektarReader::RequestData: rank: "<< my_rank<< " :  this->dynamic_coord_index = "<< this->dynamic_coord_index<< "  this->curObj->dynamic_mesh= "<<this->curObj->dynamic_mesh <<" this->DynamicMesh= "<<this->DynamicMesh);
      if(this->dynamic_coord_index >=0 && ( (this->curObj->dynamic_mesh != this->DynamicMesh) || (this->curObj->dynamic_mesh_scale != this->DynamicMeshScale ) ) )
      {
	      vtkDebugMacro(<< "vtkNektarReader::RequestData: rank: "<< my_rank<< " : set CALC_{BOUNDARY_}GEOM_FLAG to true");
	 this->CALC_GEOM_FLAG = true;
	 this->CALC_BOUNDARY_GEOM_FLAG = true;
      }
      // if the resolutions are not the same, also need to calculate the continuum geometry
      if(this->curObj->resolution != this->ElementResolution)
        {
        this->CALC_GEOM_FLAG = true;
        }
      // if the boundary resolutions are not the same, also need to calculate the boundary geometry
      if(this->curObj->boundary_resolution != this->BoundaryResolution)
        {
        this->CALC_BOUNDARY_GEOM_FLAG = true;
        }
      }

    //this->I_HAVE_DATA = false;
    }

  vtkDebugMacro(<< "vtkNektarReader::RequestData:  ElementResolution: " << ElementResolution );

  // if I have not yet read the geometry, this should only happen once
  if(this->READ_GEOM_FLAG)
    {
    if(true == vtkNektarReader::NEED_TO_MANAGER_INIT)
      {
      vtkDebugMacro( << "vtkNektarReader::RequestData: READ_GEOM_FLAG==true, call manager_init()" );

      manager_init();       /* initialize the symbol table manager */
      vtkDebugMacro(<< "vtkNektarReader::RequestData: manager_init complete" );
      vtkNektarReader::NEED_TO_MANAGER_INIT = false;
      }
    else
      {
      vtkDebugMacro( << "vtkNektarReader::RequestData: READ_GEOM_FLAG==true, NO manager_init(), already called" );
      }
    option_set("FEstorage",1);
    option_set("iterative",1);
    option_set("GSLEVEL", min(vtkMultiProcessController::GetGlobalController()->GetNumberOfProcesses(),8));
    iparam_set("NORDER.REQ", UNSET);
    iparam_set("VIS_RES", this->ElementResolution);

    //read mesh file and header of the  data file
    this->setActive();
    vtkNew<vtkTimerLog> timer;
    vtkDebugMacro( << "vtkNektarReader::RequestData: call setup with this->USE_MESH_ONLY= "<<this->USE_MESH_ONLY );
    this->nfields = setup (&this->fl, this->master, &this->nfields, 1, this->USE_MESH_ONLY);
    if(this-> NEED_TO_BACKUP_MESH)
    {
	    Store_static_mesh(this->master[0]);
    }
    
    timer->StopTimer();
    timer_diff = timer->GetElapsedTime();

    this->READ_GEOM_FLAG = false;

    vtkDebugMacro(<<"vtkNektarReader::RequestData:: Rank: "<<my_rank
                  <<" ::Done calling setup (read mesh: "<< this->p_rea_file<<"): Setup Time: "<< timer_diff);
    } // if(this->READ_GEOM_FLAG)

  if(this->READ_BOUNDARY_GEOM_FLAG)
    {
    iparam_set("BOUNDARY_RES", this->BoundaryResolution);
    // begin: NEW Boundary call
    this->setActive();
    vtkDebugMacro(<<"vtkNektarReader::RequestData: Rank: "<< my_rank<<" Active: "
                  << this->my_patch_id << " : About to call ReadBCs, this->Ubc = "<< this->Ubc);

    this->Ubc = ReadBCs(this->fl.rea.fp, this->master[0]->fhead);
    vtkDebugMacro(<<"vtkNektarReader::RequestData: Rank: "<< my_rank
                  <<" Done with call to ReadBCs, this->Ubc = "<< this->Ubc);

    for(Bndry* B=this->Ubc; B; B = B->next)
      {
      if(B->type == 'W' || B->type == 'L' || B->type == 's' || B->type == 'G' || B->type == 'C' || B->type == 'A')
        {
        B->elmt->Surface_geofac(B);
        }
      }

    this->READ_BOUNDARY_GEOM_FLAG = false;
    this->HAVE_BOUNDARY_GEOM_FLAG = true;
    } // if(this->READ_GEOM_FLAG)

  // if we need to read the data from disk..

  if(!this->I_HAVE_DATA && !this->USE_MESH_ONLY)
    {

    //sprintf(this->p_rst_file, this->p_rst_format, this->p_rst_current);
    sprintf(this->p_rst_file, this->p_rst_format, this->requested_step);
    sprintf(this->fl.in.name,"%s", this->p_rst_file);
    this->fl.in.fp = fopen(this->fl.in.name,"r");
    if (!this->fl.in.fp)
      {
      error_msg(Restart: no dumps read from restart file);
      }

    /* read the field data */
    if(FIRST_DATA)
      {
      vtkDebugMacro(<<"vtkNektarReader::RequestData: Rank: "<< my_rank<<" Now reading data from file: "<< this->p_rst_file);
      this->setActive();
      vtkNew<vtkTimerLog> timer;
      timer->StartTimer();
      ReadCopyField(&this->fl,this->master);

      timer->StopTimer();
      timer_diff = timer->GetElapsedTime();

      vtkDebugMacro(<<"vtkNektarReader::RequestData:: Rank: "<<my_rank<<" ::Done reading data from file: "<< this->p_rst_file<<":: Read  time: "<< timer_diff);
      this->curObj->setDataFilename(this->p_rst_file);
      FIRST_DATA = false;
      }
    else
      {
      this->setActive();
      vtkNew<vtkTimerLog> timer;
      timer->StartTimer();
      ReadAppendField(&this->fl,this->master, 0);
      timer->StopTimer();
      timer_diff = timer->GetElapsedTime();
      vtkDebugMacro(<<"vtkNektarReader::RequestData: Rank: "<<my_rank<<" ::Done reading (Append) data from file: "<< this->p_rst_file<<":: Read  time: "<< timer_diff);
      this->curObj->setDataFilename(this->p_rst_file);
      }
 

    /* transform field into physical space */
    vtkNew<vtkTimerLog> timer;
    for(int i = 0; i < this->nfields; ++i)
      {
      this->master[i]->Trans(this->master[i],J_to_Q);
      }
    timer->StopTimer();
    timer_diff = timer->GetElapsedTime();
    vtkDebugMacro(<<"vtkNektarReader::RequestData: Rank: "<<my_rank<<" :: Transform field into physical space time: "<< timer_diff);
    this->I_HAVE_DATA = true;
    this->memory_step = this->requested_step;

    iparam_set("VIS_RES", this->ElementResolution);
    iparam_set("BOUNDARY_RES", this->BoundaryResolution);
    } // if(!this->I_HAVE_DATA && !this->USE_MESH_ONLY)
  else
    {
    iparam_set("VIS_RES", this->ElementResolution);
    iparam_set("BOUNDARY_RES", this->BoundaryResolution);
    }

  vtkNew<vtkTimerLog> timer;
  timer->StartTimer();
  this->updateVtuData(ugrid, boundary_ugrid); // , outputPort);
  timer->StopTimer();
  timer_diff = timer->GetElapsedTime();
  vtkDebugMacro(<<"vtkNektarReader::RequestData: Rank: "<<my_rank<<" :: updateVtuData time: "<<  timer_diff);

  if(!this->USE_MESH_ONLY)
    {
    this->SetDataFileName(this->curObj->dataFilename);
    }
  total_timer->StopTimer();
  total_timer_diff = total_timer->GetElapsedTime();

  vtkDebugMacro(<<"vtkNektarReader::RequestData: Rank: "<<my_rank<< "  outputPort: " << outputPort <<" EXIT :: Total time: "<< total_timer_diff);
  return 1;
}

void vtkNektarReader::updateVtuData(vtkUnstructuredGrid* pv_ugrid, vtkUnstructuredGrid* pv_boundary_ugrid) //, int outputPort)
{
  register int i,j,k,n,e,nelmts;
  int      qa,boundary_qa,cnt;
  int      alloc_res;
  const int    nel = this->master[0]->nel;
  int      dim = this->master[0]->fhead->dim(),ntot;
  double   *z,*w, ave;
  char     *outformat;
  double timer_diff;
  double cal_stress_tensor_timer_diff;
  double interpolate_timer_diff;

  int my_rank = vtkMultiProcessController::GetGlobalController()->GetLocalProcessId();
  int num_procs = vtkMultiProcessController::GetGlobalController()->GetNumberOfProcesses();
  vtkDebugMacro(<<"vtkNektarReader::updateVtuData: ENTER: Rank: "<<my_rank); //<<" :: outputPort: "<< outputPort);
  // if the grid in the curObj is not NULL, we may have everything we need
  if(this->curObj->ugrid)
    {
    vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": this->curObj->ugrid != NULL, see if it matches");
    // if the curObj matches the request, just shallow copy, and we're done
//     if( (curObj->resolution == this->ElementResolution) &&
//         (curObj->boundary_resolution == this->BoundaryResolution) &&
//         (curObj->use_projection  == this->UseProjection) &&
//         (this->GetPointArrayStatus("Pressure") == this->curObj->pressure) &&
//         (this->GetPointArrayStatus("Velocity") == this->curObj->velocity) &&
//         (this->GetDerivedVariableArrayStatus("Vorticity")  == this->curObj->vorticity) &&
//         (this->GetDerivedVariableArrayStatus("lambda_2") == this->curObj->lambda_2) &&
//         (this->GetDerivedVariableArrayStatus("Wall Shear Stress") == this->curObj->wss ) && 
// 	(this->HAVE_BOUNDARY_GEOM_FLAG && this->curObj->boundary_ugrid->GetNumberOfCells()>0) )
    if(this->objectMatchesRequest())  
      {
          // copy the ugrid
	  pv_ugrid->ShallowCopy(this->curObj->ugrid);

	  // check to see if we should copy the bounday_ugrid
	  if(this->GetDerivedVariableArrayStatus("Wall Shear Stress") || this->HAVE_BOUNDARY_GEOM_FLAG)
	  {
		  // by the time we get here, this should not be true, this if/print is for debugging..	  
	      if(NULL == this->curObj->boundary_ugrid)
	      {
		  vtkDebugMacro(<<"vtkNektarReader::updateVtuData: NULL == this->curObj->boundary_ugrid");
	      }
	      else 
	      {
		  vtkIdType num_cells = this->curObj->boundary_ugrid->GetNumberOfCells();
		  vtkDebugMacro(<<"vtkNektarReader::updateVtuData: NULL != this->curObj->boundary_ugrid:  num_cells = "<<num_cells);
		  
	      }
	      pv_boundary_ugrid->ShallowCopy(this->curObj->boundary_ugrid);
	  }
	  else
	  {
	      if( pv_boundary_ugrid)
	      {
		  pv_boundary_ugrid->Initialize();
	      }
	  }
	  this->displayed_step = this->requested_step;
	  vtkDebugMacro(<<"vtkNektarReader::updateVtuData: ugrid same, copy : Rank: "<<my_rank); //<<" :: outputPort: "<< outputPort);
	  if(!this->USE_MESH_ONLY)
	  {
	      this->SetDataFileName(curObj->dataFilename);
	  }
	  return;
      }
    // else if the request is for less than what is in the curObj,
    // remove unwanted, update the curObj, shallow copy and we're done.
//     else if( (this->curObj->resolution == this->ElementResolution) &&
//              (this->curObj->boundary_resolution == this->BoundaryResolution) &&
//              (this->curObj->use_projection == this->UseProjection) &&
//              (!this->GetPointArrayStatus("Pressure") ||
//               (this->GetPointArrayStatus("Pressure") && this->curObj->pressure)) &&
//              (!this->GetPointArrayStatus("Velocity") ||
//               (this->GetPointArrayStatus("Velocity") && this->curObj->velocity)) &&
//              (!this->GetDerivedVariableArrayStatus("Vorticity") ||
//               (this->GetDerivedVariableArrayStatus("Vorticity") && this->curObj->vorticity)) &&
//              (!this->GetDerivedVariableArrayStatus("lambda_2") ||
//               (this->GetDerivedVariableArrayStatus("lambda_2") && this->curObj->lambda_2)) &&
//              (!this->GetDerivedVariableArrayStatus("Wall Shear Stress") ||
//               (this->GetDerivedVariableArrayStatus("Wall Shear Stress") && this->curObj->wss)) &&
// 	     (!this->HAVE_BOUNDARY_GEOM_FLAG || 
// 	      (this->HAVE_BOUNDARY_GEOM_FLAG && this->curObj->boundary_ugrid->GetNumberOfCells()>0)) )
    else if(this->objectHasExtraData())
      {
      // Remove pressure if curObj has it, but not in request
//       if(!this->GetPointArrayStatus("Pressure") && this->curObj->pressure)
//         {
//         // Does PV already have this array?  If so, remove it.
//         if (pv_ugrid->GetPointData()->GetArray("Pressure") != NULL)
//           {
//           pv_ugrid->GetPointData()->RemoveArray("Pressure");
//           }
//         // Do I already have this array?  If so, remove it.
//         if (this->curObj->ugrid->GetPointData()->GetArray("Pressure") != NULL)
//           {
//           this->curObj->ugrid->GetPointData()->RemoveArray("Pressure");
//           }
//         this->curObj->pressure = false;
//         }
//       // Remove velocity if curObj has it, but not in request
//       if(!this->GetPointArrayStatus("Velocity") && this->curObj->velocity)
//         {
//         // Does PV already have this array?  If so, remove it.
//         if (pv_ugrid->GetPointData()->GetArray("Velocity") != NULL)
//           {
//           pv_ugrid->GetPointData()->RemoveArray("Velocity");
//           }
//         // Do I already have this array?  If so, remove it.
//         if (this->curObj->ugrid->GetPointData()->GetArray("Velocity") != NULL)
//           {
//           this->curObj->ugrid->GetPointData()->RemoveArray("Velocity");
//           }
//         this->curObj->velocity = false;
//         }

	      
       //  For each of the stored variables, if it is in the curObj, but not in the request, remove it
	      for(int vid=0; vid<this->num_vars; vid++)
	      {
		      if(!this->GetPointArrayStatus(vid) && this->curObj->vars[vid])
		      {
			      // Does PV already have this array?  If so, remove it.
			      if (pv_ugrid->GetPointData()->GetArray(this->var_names[vid]) != NULL)
			      {
				      pv_ugrid->GetPointData()->RemoveArray(this->var_names[vid]);
			      }
			      // Do I already have this array?  If so, remove it.
			      if (this->curObj->ugrid->GetPointData()->GetArray(this->var_names[vid]) != NULL)
			      {
				      this->curObj->ugrid->GetPointData()->RemoveArray(this->var_names[vid]);
			      }
			      this->curObj->vars[vid] = false;
		      }
	      }
      
    
      // Remove vorticity if curObj has it, but not in request
      if(!this->GetDerivedVariableArrayStatus("Vorticity") && this->curObj->vorticity)
        {
        // Does PV already have this array?  If so, remove it.
        if (pv_ugrid->GetPointData()->GetArray("Vorticity") != NULL)
          {
          pv_ugrid->GetPointData()->RemoveArray("Vorticity");
          }
        // Do I already have this array?  If so, remove it.
        if (this->curObj->ugrid->GetPointData()->GetArray("Vorticity") != NULL)
          {
          this->curObj->ugrid->GetPointData()->RemoveArray("Vorticity");
          }
        this->curObj->vorticity = false;
        }
      // Remove lambda_2 if curObj has it, but not in request
      if(!this->GetDerivedVariableArrayStatus("lambda_2") && this->curObj->lambda_2)
        {
        // Does PV already have this array?  If so, remove it.
        if (pv_ugrid->GetPointData()->GetArray("lambda_2") != NULL)
          {
          pv_ugrid->GetPointData()->RemoveArray("lambda_2");
          }
        // Do I already have this array?  If so, remove it.
        if (this->curObj->ugrid->GetPointData()->GetArray("lambda_2") != NULL)
          {
          this->curObj->ugrid->GetPointData()->RemoveArray("lambda_2");
          }
        this->curObj->lambda_2 = false;
        }
       // Remove stress tensor if curObj has it, but not in request
      if(!this->GetDerivedVariableArrayStatus("Stress Tensor") && this->curObj->stress_tensor)
        {
        // Does PV already have this array?  If so, remove it.
        if (pv_ugrid->GetPointData()->GetArray("Stress Tensor") != NULL)
          {
          pv_ugrid->GetPointData()->RemoveArray("Stress Tensor");
          }
	if (pv_ugrid->GetPointData()->GetArray("Stress Tensor Magnitude") != NULL)
          {
          pv_ugrid->GetPointData()->RemoveArray("Stress Tensor Magnitude");
          }
        // Do I already have this array?  If so, remove it.
        if (this->curObj->ugrid->GetPointData()->GetArray("Stress Tensor") != NULL)
          {
          this->curObj->ugrid->GetPointData()->RemoveArray("Stress Tensor");
          }
	if (this->curObj->ugrid->GetPointData()->GetArray("Stress Tensor Magnitude") != NULL)
          {
          this->curObj->ugrid->GetPointData()->RemoveArray("Stress Tensor Magnitude");
          }
        this->curObj->stress_tensor = false;
        }

      pv_ugrid->ShallowCopy(this->curObj->ugrid);

      if(!this->GetDerivedVariableArrayStatus("Wall Shear Stress") && this->curObj->wss)
        {
        //pv_wss_ugrid->ShallowCopy(this->curObj->ugrid);
	    vtkDebugMacro(<<"vtkNektarReader::updateVtuData:: Rank= "<<my_rank<< " : !this->GetDerivedVariableArrayStatus(Wall Shear Stress) && this->curObj->wss");
	    if(this->curObj->boundary_ugrid == NULL)
	    {
		vtkDebugMacro(<<"vtkNektarReader::updateVtuData::this->curObj->boundary_ugrid == NULL");
	    }
	    else
	    {
		vtkDebugMacro(<<"vtkNektarReader::updateVtuData::this->curObj->boundary_ugrid != NULL");
	    }
        this->curObj->boundary_ugrid = vtkUnstructuredGrid::New();
        pv_boundary_ugrid->ShallowCopy(this->curObj->boundary_ugrid);
        this->curObj->wss = false;
        }
      else
        {
        pv_boundary_ugrid->ShallowCopy(this->curObj->boundary_ugrid);
        }
      this->displayed_step = this->requested_step;
      if(!this->USE_MESH_ONLY)
        {
        this->SetDataFileName(curObj->dataFilename);
        }
      return;
      }
    }

  // if the boundary_ugrid is NULL, we'll need to calculate WSS (before interpolating).  this is likely a new time step.
  bool NEED_TO_CALC_WSS = false;
  if(this->curObj->boundary_ugrid == NULL || ( this->boundary_mem_step != this->requested_step))
    {
    vtkDebugMacro(<<"+++++++++++++++++++++++++++  need to recalc WSS");
    NEED_TO_CALC_WSS = true;
    }
  else
    {
    vtkDebugMacro(<<"---------------------------  NO need to recalc WSS");
    }

  // otherwise the grid in the curObj is NULL, and/or the resolution has changed,
  // and/or we need more data than is in curObj, we need to do everything

  vtkDebugMacro(<<"vtkNektarReader::updateVtuData:: can't use existing ugrid, boundary_ugrid, or both, Rank = "<<my_rank);
  // <<" :: outputPort: "<< outputPort<<": this->master[0]->nel = "<< this->master[0]->nel);

  int Nvert_total = 0;
  int Nelements_total = 0;
  int vort_index = 0;
  int lambda_index = 0;
  int stress_tensor_index = 0;

  int num_total_boundary_elements = 0;
  int num_total_boundary_verts = 0;

  //int vert_ID_array_length = 0;
  int Nel = this->master[0]->nel;
  vtkSmartPointer<vtkPoints> points;
  vtkSmartPointer<vtkPoints> boundary_points;

  if(this->master[0]->fhead->dim() == 3)
    {
	    // temp_array will be used to store values for one element.  needs to be big enough for all fields. 
	    //  Later will be use for calculating vorticity and lambda_2, which will only use the first 4 values.
    double ***num;
    double *temp_array[this->nfields];

    qa = iparam("VIS_RES");
    boundary_qa = iparam("BOUNDARY_RES");
    vtkDebugMacro(<<"vtkNektarReader::updateVtuData:: rank = "<<my_rank<<" QGmax = "<<QGmax<<", qa = "<< qa << ", boundary_qa = "<< boundary_qa);
    alloc_res = (QGmax > qa) ? QGmax : qa;
    alloc_res = (alloc_res > boundary_qa) ? alloc_res : boundary_qa;
    vtkDebugMacro(<<"vtkNektarReader::updateVtuData:: rank = "<<my_rank<<" alloc_res = "<<alloc_res<<", alloc^3-1= "<<(alloc_res*alloc_res*alloc_res-1));
    Element* F  = this->master[0]->flist[0];
    vtkDebugMacro(<<"vtkNektarReader::updateVtuData:: rank = "<<my_rank<<" : F->qa = "<< F->qa);

    if(!this->USE_MESH_ONLY)
      {
        for(i=0; i<this->nfields; i++)
        {
        temp_array[i]= dvector(0,alloc_res*alloc_res*alloc_res-1);
        }
      vtkDebugMacro(<<"vtkNektarReader::updateVtuData:: rank = "<<my_rank<<" : DONE allocate temp_array");
      }

/* Paris: Adding support to other element types */
#if 0
    // calculate the number of elements and vertices in the continuum mesh
    int verts_per_element = qa*(qa+1)*(qa+2)/6;
    int sub_els_per_element = (qa-1)*(qa-1)*(qa-1);
    
    for(k = 0,n=i=0; k < Nel; ++k)
    {
      i += verts_per_element;
      n += sub_els_per_element;
      
    //F  = this->master[0]->flist[k];
    //if(Check_range_sub_cyl(F))
    //  {
      //if(Check_range(F)){
      //qa = F->qa;
      //i += qa*(qa+1)*(qa+2)/6;
      //n += (qa-1)*(qa-1)*(qa-1);
    //  }
      
    }
#endif

   for(k = 0,n=i=0; k < Nel; ++k){
      switch(F->identify()){
        case Nek_Tet:
          i += qa*(qa+1)*(qa+2)/6;
          n += (qa-1)*(qa-1)*(qa-1);
        break;
        case Nek_Hex:
          i += qa*qa*qa;
          n += (qa-1)*(qa-1)*(qa-1);
        break;
        case Nek_Prism:
          i += qa*qa*(qa+1)/2;
          n += (qa-1)*(qa-1)*(qa-1);
        break;
        default:
          fprintf(stderr,"Element is not setup in updateVtuData\n");
          exit(-1);
        break;
      }
    }

    Nvert_total = i;
    Nelements_total = n;
    vtkDebugMacro(<<"updateVtuData: rank = "<<my_rank<<" :Nvert_total= "<<Nvert_total<<", Nelements_total= "<<Nelements_total<<", Nel= "<< Nel);


/* Paris: Adding support to other element types */
#if 0
    // Calculate total number of vertices and elements for the Boundary mesh
    for(Bndry* B=this->Ubc; B; B = B->next)
      {
      if(B->type == 'W' || B->type == 'L' || B->type == 's' || B->type == 'G')
        {
        // boundary_qa is set to BOUNDARY_RES
        num_total_boundary_elements += (boundary_qa-1)*(boundary_qa-1);
        num_total_boundary_verts += (boundary_qa+1)*boundary_qa/2;
        }
      }
#else
    for(Bndry* B=this->Ubc; B; B = B->next)
      {
      if(B->type == 'W' || B->type == 'L' || B->type == 's' || B->type == 'G' || B->type == 'C' || B->type == 'A')
        {
        // boundary_qa is set to BOUNDARY_RES
         if(F->identify() == Nek_Tet){
              num_total_boundary_elements += (boundary_qa-1)*(boundary_qa-1);
              num_total_boundary_verts += (boundary_qa+1)*boundary_qa/2;
         }
	 else if(F->identify() == Nek_Hex){
              num_total_boundary_elements += (boundary_qa-1)*(boundary_qa-1);
              num_total_boundary_verts += boundary_qa*boundary_qa;
         }
         else if(F->identify() == Nek_Prism){
              num_total_boundary_elements += (boundary_qa-1)*(boundary_qa-1);
              num_total_boundary_verts += (boundary_qa+1)*boundary_qa/2;
         }
      }
    }
#endif
#if 0
    /* Paris: add cases here ?? */
    int element_vert_cnt;
    F  = this->Ubc->elmt;
    if(F->Nfverts(B->face) == 3) /* Tets and Prisms */
        element_vert_cnt = (boundary_qa-1)*(boundary_qa-1)*3;
    else 
        element_vert_cnt = (boundary_qa-1)*(boundary_qa-1); 
#endif
    /* Paris: I think not...just allocate enough memory */
 
    int element_vert_cnt = (boundary_qa-1)*(boundary_qa-1)*3;   
    int boundary_index[element_vert_cnt];


    // Generate the connectivity for one Boundary element
    generateBoundaryConnectivity(boundary_index, boundary_qa);

    // declare and allocate memory for needed scalars/vectors
    // this will need to be made into variable length arrays to accommodate
    // arbitrary numbers of variables in the next gen code

    // these should each probably be wrapped in an if statement, checking whether user requested the data

    vtkDebugMacro(<<"updateVtuData:: rank = "<<my_rank<<" : num_total_boundary_elements= "<<num_total_boundary_elements<<", num_total_boundary_verts= "<<num_total_boundary_verts);
    vtkDebugMacro(<<"updateVtuData:: rank = "<<my_rank<<" : boundary_qa= "<<boundary_qa<<", (boundary_qa+1)*boundary_qa/2= "<<((boundary_qa+1)*boundary_qa/2));

    
    vtkSmartPointer<vtkFloatArray> vorticity;
    vtkSmartPointer<vtkFloatArray> lambda_2;
    vtkSmartPointer<vtkFloatArray> stress_tensor;
    vtkSmartPointer<vtkFloatArray> stress_tensor_mag;
    vtkSmartPointer<vtkFloatArray> wall_shear_stress;

    if(!this->USE_MESH_ONLY)
      {
      if(this->GetDerivedVariableArrayStatus("Vorticity"))
        {
        vorticity = vtkSmartPointer<vtkFloatArray>::New();
        vorticity->SetNumberOfComponents(3);
        vorticity->SetNumberOfTuples(Nvert_total);
        vorticity->SetName("Vorticity");
        }

      if(this->GetDerivedVariableArrayStatus("Stress Tensor"))
        {
        stress_tensor = vtkSmartPointer<vtkFloatArray>::New();
        stress_tensor->SetNumberOfComponents(9);
        stress_tensor->SetNumberOfTuples(Nvert_total);
        stress_tensor->SetName("Stress Tensor");

	stress_tensor_mag = vtkSmartPointer<vtkFloatArray>::New();
        stress_tensor_mag->SetNumberOfComponents(1);
        stress_tensor_mag->SetNumberOfValues(Nvert_total);
        stress_tensor_mag->SetName("Stress Tensor Magnitude");
        }

      if(this->GetDerivedVariableArrayStatus("lambda_2"))
        {
        lambda_2 = vtkSmartPointer<vtkFloatArray>::New();
        lambda_2->SetNumberOfComponents(1);
        lambda_2->SetNumberOfValues(Nvert_total);
        lambda_2->SetName("lambda_2");
        }

      if(this->GetDerivedVariableArrayStatus("Wall Shear Stress"))
        {
        wall_shear_stress = vtkSmartPointer<vtkFloatArray>::New();
        wall_shear_stress->SetNumberOfComponents(3);
        wall_shear_stress->SetNumberOfValues(num_total_boundary_verts);
        wall_shear_stress->SetName("Wall Shear Stress");
        }
      } // if(!this->USE_MESH_ONLY)

    // if we need to calculate the geometry (first time, or it has changed)
    if (this->CALC_GEOM_FLAG)
      {
      // first for the full (continuum) mesh

      vtkNew<vtkTimerLog> timer;
      timer->StartTimer();
      if(this->UGrid)
        {
        this->UGrid->Delete();
        }
      this->UGrid = vtkUnstructuredGrid::New();
      this->UGrid->Allocate(Nelements_total);

      points = vtkSmartPointer<vtkPoints>::New();
      points->SetNumberOfPoints(Nvert_total);

      vtkDebugMacro(<<"updateVtuData : rank = "<<my_rank<<": Nelements_total = "<<Nelements_total<<" Nvert_total = "<< Nvert_total);

      /* fill XYZ  arrays */
      //index = 0;
      interpolateAndCopyContinuumPoints(alloc_res, qa, points);

      //vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<" number of points added:"<< index);
      timer->StopTimer();
      timer_diff = timer->GetElapsedTime();
      vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": time to copy/convert xyz and uvw: "<< timer_diff);
      } // if (this->CALC_GEOM_FLAG)

    if (this->CALC_BOUNDARY_GEOM_FLAG)
      {
      // now for the Boundary mesh (if needed)

      //if(this->GetDerivedVariableArrayStatus("Wall Shear Stress"))
        {
        vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": *** GetDerivedVariableArrayStatus(Wall Shear Stress)");

        if(this->Boundary_UGrid)
          {
          this->Boundary_UGrid->Delete();
          }
        this->Boundary_UGrid = vtkUnstructuredGrid::New();

        this->Boundary_UGrid->Allocate(num_total_boundary_elements);

        boundary_points = vtkSmartPointer<vtkPoints>::New();
        boundary_points->SetNumberOfPoints(num_total_boundary_verts);

        vtkDebugMacro(<<"updateVtuData:: my_rank= " << my_rank<<" : post SetNumberOfPoints("<<num_total_boundary_verts<<")");

	interpolateAndCopyBoundaryPoints(alloc_res, boundary_qa, boundary_points);
        vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": post interpolateAndCopyBoundaryPoints()");
	
        } // if(this->GetDerivedVariableArrayStatus("Wall Shear Stress"))
      } // if (this->CALC_BOUNDARY_GEOM_FLAG)

    if(!this->USE_MESH_ONLY)
      {
      vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": call interpolateAndCopyContinuumData()");
      interpolateAndCopyContinuumData(pv_ugrid, temp_array, qa, Nvert_total);

      vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": finished with interpolateAndCopyContinuumData (pressure and velocity) , see if we want vorticity and/or lambda_2");

      if(this->GetDerivedVariableArrayStatus("Vorticity") || this->GetDerivedVariableArrayStatus("lambda_2"))
        {
        vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": Yes, we want vorticity and/or lambda_2");
        double *vel_array[4];
        vtkNew<vtkTimerLog> timer;
        timer->StartTimer();
        for (k=0; k < 3; k++)
          {
          vel_array[k]= dvector(0, master[this->velocity_index + k]->htot*sizeof(double));
          memcpy(vel_array[k], master[this->velocity_index +k]->base_h, master[this->velocity_index +k]->htot*sizeof(double));
          }

	vel_array[k]= dvector(0, master[this->pressure_index]->htot*sizeof(double));
	memcpy(vel_array[k], master[this->pressure_index]->base_h, master[this->pressure_index]->htot*sizeof(double));
	  
        timer->StopTimer();
        timer_diff = timer->GetElapsedTime();
        vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": time to backup master: "<< timer_diff);
        vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": call Calc_Vort()");
        timer->StartTimer();
        /// we may be able to do this only once and store the values... look at this later

        if(this->UseProjection)
          option_set("PROJECT", true);
        else
          option_set("PROJECT", false);
        Calc_Vort(&this->fl, this->master, this->nfields, this->velocity_index, this->pressure_index, 0);
        timer->StopTimer();
        timer_diff = timer->GetElapsedTime();
        vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": Calc_Vort() complete, this->nfields: "<<this->nfields<<" :: time: "<< timer_diff );
        vort_index = 0;
	lambda_index = 0;
        for(k = 0; k < Nel; ++k)
          {
          F  = this->master[0]->flist[k];
          
            // 	Here we are only interpolating for vorticity and lambda_2, other values were interpolated and copied already.
		    
		    // ** jai--- where else does ntot come from?  would be better to first check if value wanted, then interpolate, then Set{Value|Tuple3}()
		    // this would mean checking ArrayStatus for each element, rather that each vertex. and only interpolating quantities that will be used.
		    
            //for(n = 0; n < this->nfields; ++n)

		    //	    verts_per_element	    
//             for(n = 0; n < 4; ++n)		    
//               {
//               F = this->master[n]->flist[k];
//               ntot = Interp_symmpts(F,qa,F->h_3d[0][0],temp_array[n],'p');
//               }

//             for(i = 0; i < ntot; ++i)
//               {
//               if(this->GetDerivedVariableArrayStatus("Vorticity"))
//                 {
//                 vorticity->SetTuple3(
//                   index,
//                   temp_array[0][i],
//                   temp_array[1][i],
//                   temp_array[2][i]);
//                 }
//               if(this->GetDerivedVariableArrayStatus("lambda_2"))
//                 {
//                 lambda_2->SetValue(index, temp_array[3][i]);
//                 }

//               index ++;
//               } // for(i = 0; i < ntot; ++i)

	  
            if(this->GetDerivedVariableArrayStatus("Vorticity"))
	    {
		    // interpolate vorticity values for this element (k)
		    for(n = 0; n < 3; ++n)		    
		    {
			    F = this->master[this->velocity_index+n]->flist[k];
			    ntot = Interp_symmpts(F,qa,F->h_3d[0][0],temp_array[n],'p');
		    }
		    // for each vertex in this element
		    for(i = 0; i < ntot; ++i)
		    {
			    vorticity->SetTuple3(
				    vort_index,
				    temp_array[0][i],
				    temp_array[1][i],
				    temp_array[2][i]);
			    vort_index++;
		    }
	    }// if(this->GetDerivedVariableArrayStatus("Vorticity"))

	    if(this->GetDerivedVariableArrayStatus("lambda_2"))
	    {
		    // interpolate lambda_2 values for this element (k)
		    n=3;
		    F = this->master[this->pressure_index]->flist[k];
		    ntot = Interp_symmpts(F,qa,F->h_3d[0][0],temp_array[n],'p');
		    
		    // for each vertex in this element
		    for(i = 0; i < ntot; ++i)
		    {
			    lambda_2->SetValue(lambda_index, temp_array[3][i]);
			    lambda_index++;
		    }
	    }// if(this->GetDerivedVariableArrayStatus("lambda_2"))

	    

	    
          }// for(k = 0; k < Nel; ++k)

        if(this->GetDerivedVariableArrayStatus("Vorticity"))
          {
          this->UGrid->GetPointData()->AddArray(vorticity);
	  //vorticity->Delete();
          }
        else
          {
          // Does PV already have this array?  If so, remove it.
          if (pv_ugrid->GetPointData()->GetArray("Vorticity") != NULL)
            {
            pv_ugrid->GetPointData()->RemoveArray("Vorticity");
            }
          // Do I already have this array?  If so, remove it.
          if (this->UGrid->GetPointData()->GetArray("Vorticity") != NULL)
            {
            this->UGrid->GetPointData()->RemoveArray("Vorticity");
            }
          }

        if(this->GetDerivedVariableArrayStatus("lambda_2"))
          {
          this->UGrid->GetPointData()->AddArray(lambda_2);
	  //lambda_2->Delete();
          }
        else  // user does not want this variable
          {
          // Does PV already have this array?  If so, remove it.
          if (pv_ugrid->GetPointData()->GetArray("lambda_2") != NULL)
            {
            pv_ugrid->GetPointData()->RemoveArray("lambda_2");
            }
          // Do I already have this array?  If so, remove it.
          if (this->UGrid->GetPointData()->GetArray("lambda_2") != NULL)
            {
            this->UGrid->GetPointData()->RemoveArray("lambda_2");
            }
          }

        timer->StartTimer();
        for (k=0; k < 3; k++)
          {
          vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": Now memcopy vel_array["<<k<<"] back to master["<<k<<"]");
          memcpy(master[this->velocity_index +k]->base_h, vel_array[k], master[this->velocity_index +k]->htot*sizeof(double));
          free(vel_array[k]);
          vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": vel_array["<<k<<"] has been freed");
          }

	memcpy(master[this->pressure_index]->base_h, vel_array[k], master[this->pressure_index]->htot*sizeof(double));
	free(vel_array[k]);
	vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": vel_array["<<k<<"] has been freed");
	  
        timer->StopTimer();
        timer_diff = timer->GetElapsedTime();
        vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": time to restore master: "<<timer_diff);

        } // if(this->GetDerivedVariableArrayStatus("Vorticity")) || this->GetDerivedVariableArrayStatus("lambda_2")
      else  // user does not want either of these variables
        {
        //if(this->GetDerivedVariableArrayStatus("Vorticity"))
        // Does PV already have this array?  If so, remove it.
        if (pv_ugrid->GetPointData()->GetArray("Vorticity") != NULL)
          {
          pv_ugrid->GetPointData()->RemoveArray("Vorticity");
          }
        // Do I already have this array?  If so, remove it.
        if (this->UGrid->GetPointData()->GetArray("Vorticity") != NULL)
          {
          this->UGrid->GetPointData()->RemoveArray("Vorticity");
          }
        // Does PV already have this array?  If so, remove it.
        if (pv_ugrid->GetPointData()->GetArray("lambda_2") != NULL)
          {
          pv_ugrid->GetPointData()->RemoveArray("lambda_2");
          }
        // Do I already have this array?  If so, remove it.
        if (this->UGrid->GetPointData()->GetArray("lambda_2") != NULL)
          {
          this->UGrid->GetPointData()->RemoveArray("lambda_2");
          }
        }

      
      if(this->GetDerivedVariableArrayStatus("Stress Tensor"))
        {
        vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": Yes, we want Stress Tensor");
        int index;
        double *backup_array[6];
	double value, FrobeniusNorm = 0.0;
        vtkNew<vtkTimerLog> timer;
        timer->StartTimer();
        for (k=0; k < 3; k++)
          {
          // allocate memory and store the displacement values
          backup_array[k]= dvector(0, master[this->sm_displacement_index + k]->htot*sizeof(double));
          memcpy(backup_array[k], master[this->sm_displacement_index +k]->base_h, master[this->sm_displacement_index +k]->htot*sizeof(double));
	  // allocate memory and store the sm velocity values (just using the storage for the values, not use in stress tensor calculation)
	  backup_array[k+3]= dvector(0, master[this->sm_velocity_index + k]->htot*sizeof(double));
          memcpy(backup_array[k+3], master[this->sm_velocity_index +k]->base_h, master[this->sm_velocity_index +k]->htot*sizeof(double));
          }

	  
        timer->StopTimer();
        timer_diff = timer->GetElapsedTime();
        vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": time to backup master (for stress tensor): "<< timer_diff);
        vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": call Calc_StressTensor(), this->sm_displacement_index: "<< this->sm_displacement_index<<"   this->sm_velocity_index: " <<  this->sm_velocity_index );
        timer->StartTimer();
        /// we may be able to do this only once and store the values... look at this later

 //        if(this->UseProjection)
//           option_set("PROJECT", true);
//         else
//           option_set("PROJECT", false);

		
        Calc_StressTensor(this->master, this->sm_displacement_index, this->sm_velocity_index);
        timer->StopTimer();
        cal_stress_tensor_timer_diff = timer->GetElapsedTime();
        vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": Calc_StressTensor() complete :: time: "<< timer_diff );

        stress_tensor_index = 0;
	
	timer->StartTimer();
        for(k = 0; k < Nel; ++k)
	{
              F  = this->master[0]->flist[k];
	  
	  
	      // interpolate stress tensor values for this element (k)
	      for(n = 0; n < 3; ++n)		    
	      {
		      F = this->master[this->sm_displacement_index+n]->flist[k];
		      ntot = Interp_symmpts(F,qa,F->h_3d[0][0],temp_array[n],'p');
		      F = this->master[this->sm_velocity_index+n]->flist[k];
		      ntot = Interp_symmpts(F,qa,F->h_3d[0][0],temp_array[n+3],'p');
		      //vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": ntot = " << ntot);
	      }
	      // for each vertex in this element
	      for(i = 0; i < ntot; ++i)
	      {
		      stress_tensor->SetTuple9(
			      stress_tensor_index,
			      temp_array[0][i],
			      temp_array[3][i],
			      temp_array[5][i],
			      temp_array[3][i],
			      temp_array[1][i],
			      temp_array[4][i],
			      temp_array[5][i],
			      temp_array[4][i],
			      temp_array[2][i]);
		      
		      for(index = 0; index < 9; ++index)
		      {
			      value = temp_array[index][i];
			      FrobeniusNorm += value*value;		  
		      }
		      FrobeniusNorm = sqrt(FrobeniusNorm);

		      stress_tensor_mag->SetValue(
			      stress_tensor_index,
			      FrobeniusNorm);
		      
		      stress_tensor_index++;
	      }
	    
	}// for(k = 0; k < Nel; ++k)

	timer->StopTimer();
	interpolate_timer_diff = timer->GetElapsedTime();

//	if(0 == my_rank)
	if(0)
	{
	    char timing_filename[256];
	    sprintf(timing_filename, "/home/insley/BROWN/FSI_%04d_procs_res_%02d.text", num_procs, this->ElementResolution);
	    FILE* timing_fd = fopen(timing_filename, "w");
	    fprintf(timing_fd, "num procs: %d \tCalc_StressTensor(): time: %f\tinterp_time: %f\ttotal: %f\n", num_procs, cal_stress_tensor_timer_diff, interpolate_timer_diff, (cal_stress_tensor_timer_diff + interpolate_timer_diff));
	    fclose(timing_fd);
	}

	//this->UGrid->GetPointData()->AddArray(stress_tensor);
	this->UGrid->GetPointData()->SetTensors(stress_tensor);
	//stress_tensor->Delete();
	this->UGrid->GetPointData()->AddArray(stress_tensor_mag);
	//stress_tensor_mag->Delete();
	vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<":  stress_tensor_index = " <<  stress_tensor_index);
	
	timer->StartTimer();
	for (k=0; k < 3; k++)
	{
          vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": Now memcopy backup_array["<<k<<"] back to master["<<this->sm_displacement_index+k<<"]");
          memcpy(master[this->sm_displacement_index +k]->base_h, backup_array[k], master[this->sm_displacement_index +k]->htot*sizeof(double));
          free(backup_array[k]);
          vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": backup_array["<<k<<"] has been freed");
	  vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": Now memcopy backup_array["<<k+3<<"] back to master["<<this->sm_velocity_index +k<<"]");
          memcpy(master[this->sm_velocity_index +k]->base_h, backup_array[k+3], master[this->sm_velocity_index +k]->htot*sizeof(double));
          free(backup_array[k+3]);
          vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": backup_array["<<k+3<<"] has been freed");
	}
	
        timer->StopTimer();
        timer_diff = timer->GetElapsedTime();
        vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": time to restore master (for stress tensor): "<<timer_diff);

        } // if(this->GetDerivedVariableArrayStatus("Stress Tensor")) 
      else  // user does not want Stress Tensor
        {
        // Does PV already have this array?  If so, remove it.
        if (pv_ugrid->GetPointData()->GetArray("Stress Tensor") != NULL)
          {
          pv_ugrid->GetPointData()->RemoveArray("Stress Tensor");
          }
        // Do I already have this array?  If so, remove it.
        if (this->UGrid->GetPointData()->GetArray("Stress Tensor") != NULL)
          {
          this->UGrid->GetPointData()->RemoveArray("Stress Tensor");
          }
        } // else (user does not want stress tensor)

      

      } // if(!this->USE_MESH_ONLY)

    // now see if they want the Wall Shear Stress (WSS)
    vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<"Should we calc WSS: "<< this->GetDerivedVariableArrayStatus("Wall Shear Stress"));

    if(this->GetDerivedVariableArrayStatus("Wall Shear Stress") || this->USE_MESH_ONLY || this->HAVE_BOUNDARY_GEOM_FLAG)
      {
      vtkIdType num_boundary_cells = this->Boundary_UGrid->GetNumberOfCells();
      if(num_boundary_cells == 0)
        {
        addCellsToBoundaryMesh(boundary_index, boundary_qa);
        this->Boundary_UGrid->SetPoints(boundary_points);
        }

      if(!this->USE_MESH_ONLY)
        {
		int some_number;	
        if(this->Ubc == NULL)
          vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": Want to call Calc_WSS(, but (this->Ubc == NULL");
        vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": Lets call Calc_WSS()");
        // if WSS_all_vals is NULL, then we need to allocate it
        if(this->WSS_all_vals == NULL)
          {
          // figure out how many total values there are
          some_number = get_number_of_boundary_vertices(&this->fl, this->master, this->Ubc, 0);

	  //vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": WSS_all_vals == NULL, so allocate for some_number= "<< some_number );
	   
          // allocate enough memory for them
          this->WSS_all_vals= new double*[3];
          for(int w=0; w<3; w++)
            {
            this->WSS_all_vals[w] = (double*) malloc(some_number *sizeof(double));
	    bzero(this->WSS_all_vals[w], some_number *sizeof(double));
            }
          }
        if(NEED_TO_CALC_WSS)
          {
          Calc_WSS(&this->fl, this->master, this->Ubc, 0, this->WSS_all_vals);
          this->boundary_mem_step = this->requested_step;
	  
          vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": Done calling Calc_WSS()");
          }

        vtkDebugMacro(<< "~~~~~updateVtuData: my_rank= " << my_rank<<":  this->Boundary_UGrid->GetNumberOfCells() = "<< num_boundary_cells);

        interpolateAndCopyBoundaryData(alloc_res, num_total_boundary_verts , boundary_qa);
        } // if(!this->USE_MESH_ONLY)

      }// if(this->GetDerivedVariableArrayStatus("Wall Shear Stress"))

    vtkNew<vtkTimerLog> timer;
    timer->StartTimer();
    if (this->CALC_GEOM_FLAG)
      {
      /* continuum numbering array */
      num = dtarray(0,alloc_res-1,0,alloc_res-1,0,alloc_res-1);
      /* Paris: Adding support to other element types */
      switch(F->identify()){
        case Nek_Tet:
          for(cnt = 1, k = 0; k < qa; ++k)
          {
            for(j = 0; j < qa-k; ++j)
            {
              for(i = 0; i < qa-k-j; ++i, ++cnt)
              {
                num[i][j][k] = cnt;
              }
            }
          }
        break;
        case Nek_Hex:
          for(cnt = 1, k = 0; k < qa; ++k)
          {
            for(j = 0; j < qa; ++j)
            {
              for(i = 0; i < qa; ++i, ++cnt)
              {
                num[i][j][k] = cnt;
              }
            }
          }
        break;
          case Nek_Prism:
            for(cnt = 1, k = 0; k < qa; ++k)
            {
              for(j = 0; j < qa; ++j)
              {
                for(i = 0; i < qa-k; ++i, ++cnt)
                {
                  num[i][j][k] = cnt;
                }
              }
           }
        break;
        default:
          fprintf(stderr,"Element is not setup in updateVtuData\n");
          exit(-1);
        break;
      }
#if 0
      for(cnt = 1, k = 0; k < qa; ++k)
        {
        for(j = 0; j < qa-k; ++j)
          {
          for(i = 0; i < qa-k-j; ++i, ++cnt)
            {
            num[i][j][k] = cnt;
            }
          }
        }
#endif
      vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<":: cnt = "<<cnt<<" : Nel = "<<Nel<<" : cnt*Nel = "<<(cnt*Nel)<< " : num= "<< num << " : num[0][0][0]= "<< num[0][0][0]);

      gsync();
      addCellsToContinuumMesh(qa, num);

      vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": vort_index = "<< vort_index << "  lambda_index = " << lambda_index);
      this->UGrid->SetPoints(points);
      free_dtarray(num,0,0,0);
      }//end "if (this->CALC_GEOM_FLAG)"

    timer->StopTimer();
    timer_diff = timer->GetElapsedTime();
    vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": time of CALC_GEOM (the mesh): "<< timer_diff);

    timer->StartTimer();
    vtkNew<vtkCleanUnstructuredGrid> clean;
#ifdef VTK5
    clean->SetInput(this->UGrid);
#else
    vtkNew<vtkUnstructuredGrid> tmpGrid;
    tmpGrid->ShallowCopy(this->UGrid);
    clean->SetInputData(tmpGrid.GetPointer());
#endif
    clean->Update();
    timer->StopTimer();
    timer_diff = timer->GetElapsedTime();
    vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": time to clean the grid: "<< timer_diff);

    timer->StartTimer();
    pv_ugrid->ShallowCopy(clean->GetOutput());

    vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<":  completed ShallowCopy to pv_ugrid\n");
    if(this->curObj->ugrid)
      {
      this->curObj->ugrid->Delete();
      }
    this->curObj->ugrid = vtkUnstructuredGrid::New();
    this->curObj->ugrid->ShallowCopy(clean->GetOutput());

    if(this->curObj->boundary_ugrid)
      {
      vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": Delete this->curObj->boundary_ugrid it is not null\n");
      this->curObj->boundary_ugrid->Delete();
      }

    this->curObj->boundary_ugrid = vtkUnstructuredGrid::New();
    vtkNew<vtkCleanUnstructuredGrid> boundary_clean;

    vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": created new this->curObj->boundary_ugrid, and boundary_clean\n");

    if(this->GetDerivedVariableArrayStatus("Wall Shear Stress") || this->USE_MESH_ONLY || this->HAVE_BOUNDARY_GEOM_FLAG)
      {
#ifdef VTK5
      boundary_clean->SetInput(this->Boundary_UGrid);
#else
      vtkNew<vtkUnstructuredGrid> tmpGrid;
      tmpGrid->ShallowCopy(this->Boundary_UGrid);
      boundary_clean->SetInputData(tmpGrid.GetPointer());
#endif
      boundary_clean->Update();
      vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": they want WSS (or boundary mesh), do shallow copy\n");
      
      this->curObj->boundary_ugrid->ShallowCopy(boundary_clean->GetOutput());
      pv_boundary_ugrid->ShallowCopy(boundary_clean->GetOutput());
      }
    else
      {
      vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": they DO NOT want WSS, do nothing\n");
      pv_boundary_ugrid->ShallowCopy(this->curObj->boundary_ugrid );
      }

    //fprintf(stderr, "updateVtuData:: Rank: %d:  completed ShallowCopy to curObj->ugrid\n");
    timer->StopTimer();
    timer_diff = timer->GetElapsedTime();
    vtkDebugMacro(<< "updateVtuData: my_rank= " << my_rank<<": time to shallow copy cleaned grid to pv_boundary_ugrid: "<<timer_diff);

    // update the current object to reflect what was requested
    this->displayed_step = this->requested_step;
    //this->curObj->pressure = this->GetPointArrayStatus("Pressure");
    // this->curObj->velocity = this->GetPointArrayStatus("Velocity");
    this->curObj->vorticity = this->GetDerivedVariableArrayStatus("Vorticity");
    this->curObj->lambda_2 = this->GetDerivedVariableArrayStatus("lambda_2");
    this->curObj->wss = this->GetDerivedVariableArrayStatus("Wall Shear Stress");
    this->curObj->stress_tensor = this->GetDerivedVariableArrayStatus("Stress Tensor");
    for(int kk=0; kk<this->num_vars; kk++)
      {
	      //this->curObj->vars[kk] = this->GetPointArrayStatus(this->var_names[kk]);
      this->curObj->vars[kk] = this->GetPointArrayStatus(kk);
      }
    this->curObj->resolution = this->ElementResolution;
    this->curObj->boundary_resolution = this->BoundaryResolution;
    this->curObj->use_projection = this->UseProjection;
    this->curObj->dynamic_mesh = this->DynamicMesh;
    this->curObj->dynamic_mesh_scale = this->DynamicMeshScale;

    if(!this->USE_MESH_ONLY)
      {
      for(i=0; i<this->nfields; i++)
        {
        free(temp_array[i]);
        }
      }
    this->CALC_GEOM_FLAG=false;
    }
  else
    {
    fprintf(stderr,"2D case is not implemented");
    }

} // vtkNektarReader::updateVtuData()

/* Paris: Adding support to other element types */
void vtkNektarReader::addCellsToContinuumMesh(int qa, double ***num)
{
  int my_rank = vtkMultiProcessController::GetGlobalController()->GetLocalProcessId();
  vtkDebugMacro(<< "addCellsToContinuumMesh(): my_rank= " << my_rank<<"  ENTER");
  vtkIdType pts[8];
  int index = 0;

  vtkDebugMacro(<< "addCellsToContinuumMesh(): my_rank= " << my_rank<< " : num= "<< num << " : num[0][0][0]= "<< num[0][0][0]);
  int n = 0;
  int nnodes;
  for(int e = 0; e < this->master[0]->nel; ++e)
    {
      //qa = F->qa;
      /* dump connectivity */
      switch(this->master[0]->flist[e]->identify())
        {
        case Nek_Tet:
          nnodes = ((qa*(qa+1)*(qa+2))/6);
          for(int k=0; k < qa-1; ++k)
            {
            for(int j = 0; j < qa-1-k; ++j)
              {
              for(int i = 0; i < qa-2-k-j; ++i)
                {

                pts[0] = n+(int) num[i][j][k] -1;
                pts[1] = n+(int) num[i+1][j][k] -1;
                pts[2] = n+(int) num[i][j+1][k] -1;
                pts[3] = n+(int) num[i][j][k+1] -1;
                this->UGrid->InsertNextCell(VTK_TETRA, 4, pts);
                index++;

                pts[0] = n+(int) num[i+1][j][k]-1;
                pts[1] = n+(int) num[i][j+1][k]-1;
                pts[2] = n+(int) num[i][j][k+1]-1;
                pts[3] = n+(int) num[i+1][j][k+1]-1;
                this->UGrid->InsertNextCell(VTK_TETRA, 4, pts);
                index++;

                pts[0] = n+(int) num[i+1][j][k+1]-1;
                pts[1] = n+(int) num[i][j][k+1]-1;
                pts[2] = n+(int) num[i][j+1][k+1]-1;
                pts[3] = n+(int) num[i][j+1][k]-1;
                this->UGrid->InsertNextCell(VTK_TETRA, 4, pts);
                index++;

                pts[0] = n+(int) num[i+1][j+1][k]-1;
                pts[1] = n+(int) num[i][j+1][k]-1;
                pts[2] = n+(int) num[i+1][j][k]-1;
                pts[3] = n+(int) num[i+1][j][k+1]-1;
                this->UGrid->InsertNextCell(VTK_TETRA, 4, pts);
                index++;

                pts[0] = n+(int) num[i+1][j+1][k]-1;
                pts[1] = n+(int) num[i][j+1][k]-1;
                pts[2] = n+(int) num[i+1][j][k+1]-1;
                pts[3] = n+(int) num[i][j+1][k+1]-1;
                this->UGrid->InsertNextCell(VTK_TETRA, 4, pts);
                index++;

                if(i < qa-3-k-j)
                  {
                  pts[0] = n+(int) num[i][j+1][k+1]-1;
                  pts[1] = n+(int) num[i+1][j+1][k+1]-1;
                  pts[2] = n+(int) num[i+1][j][k+1]-1;
                  pts[3] = n+(int) num[i+1][j+1][k]-1;
                  this->UGrid->InsertNextCell(VTK_TETRA, 4, pts);
                  index++;
                  }

                }
              pts[0] = n+(int) num[qa-2-k-j][j][k]-1;
              pts[1] = n+(int) num[qa-1-k-j][j][k]-1;
              pts[2] = n+(int) num[qa-2-k-j][j+1][k]-1;
              pts[3] = n+(int) num[qa-2-k-j][j][k+1]-1;
              this->UGrid->InsertNextCell(VTK_TETRA, 4, pts);
              index++;
              }
            }
            n += nnodes;
          break;
        case Nek_Hex:

          nnodes = qa*qa*qa;

          for(int k=0; k < qa-1; ++k)
            {
            for(int j = 0; j < qa-1; ++j)
              {
              for(int i = 0; i < qa-1; ++i)
                {
                pts[0] = n+(int) num[i][j][k] - 1;
                pts[1] = n+(int) num[i+1][j][k] - 1;
                pts[2] = n+(int) num[i+1][j+1][k] - 1;
                pts[3] = n+(int) num[i][j+1][k] -1;
                pts[4] = n+(int) num[i][j][k+1] - 1;
                pts[5] = n+(int) num[i+1][j][k+1] - 1;
                pts[6] = n+(int) num[i+1][j+1][k+1] -1;
                pts[7] = n+(int) num[i][j+1][k+1] -1;

                this->UGrid->InsertNextCell(VTK_HEXAHEDRON, 8, pts);
                index++;
              }
            }
          }
          n += nnodes;
        break;
        case Nek_Prism:

	  nnodes = qa*qa*(qa+1)/2;
          for(int k=0; k < qa-1; ++k)
            {
            for(int j = 0; j < qa-1; ++j)
              {
              for(int i = 0; i < qa-2-k; ++i)
                {

                pts[0] = n+(int) num[i][j][k] - 1;
                pts[1] = n+(int) num[i+1][j][k] - 1;
                pts[2] = n+(int) num[i][j][k+1] - 1;
                pts[3] = n+(int) num[i][j+1][k] - 1;
                pts[4] = n+(int) num[i+1][j+1][k] - 1;
                pts[5] = n+(int) num[i][j+1][k+1] - 1;

                this->UGrid->InsertNextCell(VTK_WEDGE, 6, pts);
                index++;

                pts[0] = n+(int) num[i+1][j][k] - 1;
                pts[1] = n+(int) num[i+1][j][k+1] - 1;
                pts[2] = n+(int) num[i][j][k+1] - 1;
                pts[3] = n+(int) num[i+1][j+1][k] - 1;
                pts[4] = n+(int) num[i+1][j+1][k+1] - 1;
                pts[5] = n+(int) num[i][j+1][k+1] - 1;

                this->UGrid->InsertNextCell(VTK_WEDGE, 6, pts);
                index++;

		}

              pts[0] = n+(int) num[qa-2-k][j][k]-1;
              pts[1] = n+(int) num[qa-1-k][j][k]-1;
              pts[2] = n+(int) num[qa-2-k][j][k+1]-1;
              pts[3] = n+(int) num[qa-2-k][j+1][k]-1;
              pts[4] = n+(int) num[qa-1-k][j+1][k]-1;
              pts[5] = n+(int) num[qa-2-k][j+1][k+1]-1;

              this->UGrid->InsertNextCell(VTK_WEDGE, 6, pts);
              index++;

	      }
            }
            n += nnodes;
          break;
#if 0
        case Nek_Prism:
        {
          int sum = 0;
          for (int i=qa ;i > 0; i--){
           sum += i;
          }
          nnodes = qa*sum;

          for(int k=0; k < qa-1; ++k)
            {
            for(int j = 0; j < qa-1; ++j)
              {
              for(int i = 0; i < qa-k-1; ++i)
                {
                pts[0] = n+(int) num[i][j][k] - 1;
                pts[1] = n+(int) num[i+1][j][k] - 1;
                pts[2] = n+(int) num[i][j][k+1] - 1;
                pts[3] = n+(int) num[i][j+1][k] - 1;
                this->UGrid->InsertNextCell(VTK_TETRA, 4, pts);
                index++;

                pts[0] = n+(int) num[i+1][j][k] - 1;
                pts[1] = n+(int) num[i][j][k+1] - 1;
                pts[2] = n+(int) num[i][j+1][k] - 1;
                pts[3] = n+(int) num[i+1][j+1][k] - 1;
                this->UGrid->InsertNextCell(VTK_TETRA, 4, pts);
                index++;

                pts[0] = n+(int) num[i][j][k+1] - 1;
                pts[1] = n+(int) num[i][j+1][k] - 1;
                pts[2] = n+(int) num[i+1][j+1][k] - 1;
                pts[3] = n+(int) num[i][j+1][k+1] - 1;
                this->UGrid->InsertNextCell(VTK_TETRA, 4, pts);
                index++;

                if(i < qa-k-2)
                  {
                  pts[0] = n+(int) num[i+1][j][k]-1;
                  pts[1] = n+(int) num[i+1][j][k+1]-1;
                  pts[2] = n+(int) num[i][j][k+1]-1;
                  pts[3] = n+(int) num[i+1][j+1][k]-1;
                  this->UGrid->InsertNextCell(VTK_TETRA, 4, pts);
                  index++;

                  pts[0] = n+(int) num[i+1][j][k+1]-1;
                  pts[1] = n+(int) num[i][j][k+1]-1;
                  pts[2] = n+(int) num[i+1][j+1][k]-1;
                  pts[3] = n+(int) num[i+1][j+1][k+1]-1;
                  this->UGrid->InsertNextCell(VTK_TETRA, 4, pts);
                  index++;

                  pts[0] = n+(int) num[i][j][k+1]-1;
                  pts[1] = n+(int) num[i+1][j+1][k]-1;
                  pts[2] = n+(int) num[i+1][j+1][k+1]-1;
                  pts[3] = n+(int) num[i][j+1][k+1]-1;
                  this->UGrid->InsertNextCell(VTK_TETRA, 4, pts);
                  index++;
                  }
                }
              }
          }
          n += nnodes;
        break;
        }
#endif
        default:
          fprintf(stderr,"WriteS is not set up for this element type \n");
          exit(1);
          break;
        } // switch(F->identify())
    } // for(e = 0,n=0; e <   int Nel = this->master[0]->nel; ++e)
  vtkDebugMacro(<< "addCellsToContinuumMesh(): my_rank= " << my_rank<<"  EXIT");
}// addPointsToContinuumMesh()


/* Paris: Adding support to other element types */
void vtkNektarReader::generateBoundaryConnectivity(int * boundary_index, int res)
{
  int  j, k=0;

  switch(this->master[0]->flist[0]->identify())
    {
     case Nek_Tet:
      for(int cnt_local = 0,j = 0; j < res-1; ++j)
      {
      for(int i = 0; i < res-2-j; ++i)
        {
        boundary_index[k++] = cnt_local + i + 1 - 1;
        boundary_index[k++] = cnt_local + i + 2 - 1;
        boundary_index[k++] = cnt_local + res - j + i + 1 - 1;
        boundary_index[k++] = cnt_local + res - j + i + 2 - 1;
        boundary_index[k++] = cnt_local + res - j + i + 1 - 1;
        boundary_index[k++] = cnt_local + i + 2 - 1;
        }
      boundary_index[k++] = cnt_local + res - 1 - j - 1;
      boundary_index[k++] = cnt_local + res - j - 1;
      boundary_index[k++] = cnt_local + 2*res - 2*j - 1 - 1;
      cnt_local += res - j;
      }// for(cnt_local = 0,j = 0; j < qa-1; ++j)
    break;
    case Nek_Prism: 
      for(int cnt_local = 0,j = 0; j < res-1; ++j)
      {
      for(int i = 0; i < res-2-j; ++i)
        {
        boundary_index[k++] = cnt_local + i + 1 - 1;
        boundary_index[k++] = cnt_local + i + 2 - 1;
        boundary_index[k++] = cnt_local + res - j + i + 1 - 1;
        boundary_index[k++] = cnt_local + res - j + i + 2 - 1;
        boundary_index[k++] = cnt_local + res - j + i + 1 - 1;
        boundary_index[k++] = cnt_local + i + 2 - 1;
        }
      boundary_index[k++] = cnt_local + res - 1 - j - 1;
      boundary_index[k++] = cnt_local + res - j - 1;
      boundary_index[k++] = cnt_local + 2*res - 2*j - 1 - 1;
      cnt_local += res - j;
      }// for(cnt_local = 0,j = 0; j < qa-1; ++j)
    break;
    case Nek_Hex:
      for (int cnt_local = 0, j= 0; j < res-1; ++j)
      {
      for (int i = 0; i < res-1; ++i)
        {
        boundary_index[k++] = i + cnt_local;
        boundary_index[k++] = i + cnt_local + 1;
        boundary_index[k++] = i + cnt_local + res + 1;  
        boundary_index[k++] = i + cnt_local + res;
        }
      cnt_local += res;
      }
      break;
  }

}// vtkNektarReader::generateBoundaryconnectivity()

#if 0
void vtkNektarReader::generateBoundaryConnectivity(int * boundary_index, int res)
{
  int j, k=0;

  for(int cnt_local = 0,j = 0; j < res-1; ++j)
    {
    for(int i = 0; i < res-2-j; ++i)
      {
      boundary_index[k++] = cnt_local + i + 1 - 1;
      boundary_index[k++] = cnt_local + i + 2 - 1;
      boundary_index[k++] = cnt_local + res - j + i + 1 - 1;
      boundary_index[k++] = cnt_local + res - j + i + 2 - 1;
      boundary_index[k++] = cnt_local + res - j + i + 1 - 1;
      boundary_index[k++] = cnt_local + i + 2 - 1;
      }
    boundary_index[k++] = cnt_local + res - 1 - j - 1;
    boundary_index[k++] = cnt_local + res - j - 1;
    boundary_index[k++] = cnt_local + 2*res - 2*j - 1 - 1;
    cnt_local += res - j;
    }// for(cnt_local = 0,j = 0; j < qa-1; ++j)
}// vtkNektarReader::generateBoundaryconnectivity()
#endif

void vtkNektarReader::addCellsToBoundaryMesh(int * boundary_index, int boundary_qa)
{
  int n=0;
  vtkIdType pts[4];

  Nek_Facet_Type element_type = this->master[0]->flist[0]->identify();

  for(Bndry* B=this->Ubc; B; B = B->next)
    {
    if(B->type == 'W' || B->type == 'L' || B->type == 's' || B->type == 'G' || B->type == 'C' || B->type == 'A')
      {

      int vert_index = 0;
      int cnt = 0;
      switch(element_type)
      {
        case Nek_Tet: case Nek_Prism:
          for(int j = 0; j < boundary_qa-1; ++j)
            {
            //for(i = 0; i < QGmax-2-j; ++i)
            //for(i = 0; i < alloc_res-2-j; ++i)
            for(int i = 0; i < boundary_qa-2-j; ++i)
              {
              //fprintf(out,"%d %d %d\n",cnt+i+1, cnt+i+2,cnt+boundary_qa-j+i+1);
              pts[0] = n +  boundary_index[vert_index++];
              pts[1] = n +  boundary_index[vert_index++];
              pts[2] = n +  boundary_index[vert_index++];
              this->Boundary_UGrid->InsertNextCell(VTK_TRIANGLE, 3, pts);
              //fprintf(out,"%d %d %d\n",cnt+boundary_qa-j+i+2,cnt+boundary_qa-j+i+1,cnt+i+2);
              pts[0] = n +  boundary_index[vert_index++];
              pts[1] = n +  boundary_index[vert_index++];
              pts[2] = n +  boundary_index[vert_index++];
              this->Boundary_UGrid->InsertNextCell(VTK_TRIANGLE, 3, pts);
              }
            //fprintf(out,"%d %d %d\n",cnt+boundary_qa-1-j,cnt+boundary_qa-j,cnt+2*boundary_qa-2*j-1);

            pts[0] = n +  boundary_index[vert_index++];
            pts[1] = n +  boundary_index[vert_index++];
            pts[2] = n +  boundary_index[vert_index++];
            this->Boundary_UGrid->InsertNextCell(VTK_TRIANGLE, 3, pts);
            }// for(cnt = 0,j = 0; j < alloc_res-1; ++j)

          n+= (boundary_qa+1)*boundary_qa/2;
        break;
        case Nek_Hex:
        {
          for(int j = 0; j < boundary_qa-1; ++j)
            {
            for(int i = 0; i < boundary_qa-1; ++i)
              {
              pts[0] = n +  boundary_index[vert_index++];
              pts[1] = n +  boundary_index[vert_index++];
              pts[2] = n +  boundary_index[vert_index++];
              pts[3] = n +  boundary_index[vert_index++];
              this->Boundary_UGrid->InsertNextCell(VTK_QUAD, 4, pts);
              }
            }
          n+= boundary_qa*boundary_qa;
        break; 
        }
      } // if(B->type == 'W' || B->type == 'L' || B->type == 's' || B->type == 'G')
    }// for(B=this->Ubc; B; B = B->next)
  }// addCellsToBoundaryMesh

}

void vtkNektarReader::interpolateAndCopyContinuumPoints(int alloc_res, int interp_res, vtkPoints* points)
{
  /* fill XYZ  arrays */
  Coord    X;
  int index = 0;

  int my_rank = vtkMultiProcessController::GetGlobalController()->GetLocalProcessId();
  
  // if the user requested a dynamic mesh, (and the data is available) update the mesh with the appropriate data values
  if(this->DynamicMesh && this->dynamic_coord_index >=0)
  {
	  vtkDebugMacro(<< "interpolateAndCopyContinuumCoordinates: my_rank= " << my_rank<<" Calling Set_Mesh(use dynamic) with this->dynamic_coord_index= "<< this->dynamic_coord_index);
	  Set_Mesh(this->master, this->nfields, this->dynamic_coord_index, this->DynamicMeshScale, 0);
  }
  else if(!this->DynamicMesh && this->curObj->dynamic_mesh  && this->dynamic_coord_index >=0)
  {
	  vtkDebugMacro(<< "interpolateAndCopyContinuumCoordinates: my_rank= " << my_rank<<" Calling Set_Mesh(use static) with this->dynamic_coord_index= "<< this->dynamic_coord_index);
	  Set_Mesh(this->master, this->nfields, this->dynamic_coord_index, this->DynamicMeshScale, 1);
	   vtkDebugMacro(<< "interpolateAndCopyContinuumCoordinates: my_rank= " << my_rank<<" Calling Set_Mesh(use static) complete");
  }
  
  X.x = dvector(0,alloc_res*alloc_res*alloc_res-1);
  X.y = dvector(0,alloc_res*alloc_res*alloc_res-1);
  X.z = dvector(0,alloc_res*alloc_res*alloc_res-1);

  // for each spectral element in the continuum mesh
  for(int k = 0; k < this->master[0]->nel; ++k)
    {
    Element* F  = this->master[0]->flist[k];

    //if(Check_range_sub_cyl(F))
      {
      //if(Check_range(F)){
      //qa = F->qa;

      // Get the coordinates for this element of the continuum mesh
      F->coord(&X);

      // interpolate to the requested resolution (qa == VIS_RES)
      int ntot = Interp_symmpts(F,interp_res,X.x,X.x,'p');
      ntot = Interp_symmpts(F,interp_res,X.y,X.y,'p');
      ntot = Interp_symmpts(F,interp_res,X.z,X.z,'p');

      // for every point in this spectral element
      for(int i = 0; i < ntot; ++i)
        {
        points->InsertPoint(index, X.x[i], X.y[i], X.z[i]);
        index ++;
        }// for(i = 0; i < ntot; ++i)
      }
    }//for(k = 0; k < this->master[0]->nel; ++k)

  free(X.x); free(X.y); free(X.z);
 
  vtkDebugMacro(<< "interpolateAndCopyContinuumCoordinates: my_rank= " << my_rank<<" number of points added:"<< index);
}// vtkNektarReader::interpolateContinuumCoordinates()

void vtkNektarReader::interpolateAndCopyContinuumData(vtkUnstructuredGrid* pv_ugrid, double **data_array, int interp_res, int num_verts)
{
  int my_rank = vtkMultiProcessController::GetGlobalController()->GetLocalProcessId();
  vtkDebugMacro(<< "interpolateAndCopyContinuumData: my_rank= " << my_rank<<"  ENTER");
  vtkDebugMacro(<< "interpolateAndCopyContinuumData: my_rank= " << my_rank<<"  this->num_used_scalars: " << this->num_used_scalars << "  this->num_used_vectors: "<< this->num_used_vectors);


  int index = 0;
  int cur_scalar_index=0;
  int cur_vector_index=0;
  int v_index;

//   vtkFloatArray* pressure;
//   pressure = vtkFloatArray::New();
//   pressure->SetNumberOfComponents(1);
//   pressure->SetNumberOfValues(num_verts);
//   pressure->SetName("Pressure");

  vtkFloatArray** scalars;
  scalars = (vtkFloatArray**) malloc(this->num_used_scalars * sizeof(vtkFloatArray*));

  vtkFloatArray** vectors;
  vectors = (vtkFloatArray**) malloc(this->num_used_vectors * sizeof(vtkFloatArray*));

  // allocate arrays for used scalars and vectors
  for(int jj=0; jj<this->num_vars; jj++)
    {
	if(this->GetPointArrayStatus(jj) == 1)
	{
	// if this variable is a scalar
	if(this->var_length[jj] == 1)
	{
	    vtkDebugMacro(<< "interpolateAndCopyContinuumData: my_rank= " << my_rank<<": var["<<jj<<"]: allocate scalars["<<cur_scalar_index<<"]:  name= "<< this->var_names[jj]);
	    scalars[cur_scalar_index] = vtkFloatArray::New();
	    scalars[cur_scalar_index]->SetNumberOfComponents(1);
	    scalars[cur_scalar_index]->SetNumberOfValues(num_verts);
	    scalars[cur_scalar_index]->SetName(this->var_names[jj]);
	    cur_scalar_index++;
	}
	// if this variable is a vector
	else if(this->var_length[jj] == 3)
	{
	    vtkDebugMacro(<< "interpolateAndCopyContinuumData: my_rank= " << my_rank<<": var["<<jj<<"]: allocate vectors["<<cur_vector_index<<"]:  name= "<< this->var_names[jj]);
	    vectors[cur_vector_index] = vtkFloatArray::New();
	    vectors[cur_vector_index]->SetNumberOfComponents(3);
	    vectors[cur_vector_index]->SetNumberOfTuples(num_verts);
	    vectors[cur_vector_index]->SetName(this->var_names[jj]);
	    cur_vector_index++;
	}
	}
//   vtkFloatArray** vectors;
//   vectors = vtkFloatArray::New();
//   vectors->SetNumberOfComponents(3);
//   vectors->SetNumberOfTuples(num_verts);
//   vectors->SetName("Velocity");
    }

  // for each spectral element in the continuum mesh
  for(int k = 0; k < this->master[0]->nel; ++k)
    {
    Element* F  = this->master[0]->flist[k];
      {
      // interpolate the field data on the requested resolution, storing result in data_array
      int ntot = 0;
      for(int n = 0; n < this->nfields; ++n)
      {
	  F = this->master[n]->flist[k];
	  // only interpolate fields that we need
	  if(this->use_field[n])
	  {
	      //vtkDebugMacro(<< "interpolateAndCopyContinuumData: my_rank= " << my_rank<<"  Interp_symmpts for field "<< n);
	      ntot = Interp_symmpts(F,interp_res,F->h_3d[0][0],data_array[n],'p');
	  }
      }
      
 
      
      // it may be faster to look at each variable first, and then walk through all elements then all points
      // need to look into this

      // for every point in this spectral element
      for(int i = 0; i < ntot; ++i)
      {
	  int field_index=0;
	  cur_scalar_index=0;
	  cur_vector_index=0;

	  // for each variable
	  for(v_index=0; v_index< this->num_vars; v_index++)
	  {
	      // if the variable is requested, add it to array of scalars or vectors, accordingly
	      //if(this->GetPointArrayStatus(v_index) == 1)
	      //{

	      // if this is a scalar
	      if(this->var_length[v_index] == 1)
	      {
		  // if it is requested, add it to scalars array
		  if(this->GetPointArrayStatus(v_index) == 1)
		  {
		      scalars[cur_scalar_index++]->SetValue(index, data_array[field_index][i]);
		  }
		  field_index++;
	      }
	      // if this is a vector
	      else if(this->var_length[v_index] == 3)
	      {
		  // if it is requested, add it to vectors array
		  if(this->GetPointArrayStatus(v_index) == 1)
		  {
		      vectors[cur_vector_index++]->SetTuple3(index, 
							     data_array[field_index][i],
							     data_array[field_index+1][i],
							     data_array[field_index+2][i]);
		  }
		  field_index+=3;
	      }
	      else
	      {
		  vtkDebugMacro(<< "interpolateAndCopyContinuumData: my_rank= " << my_rank<<" :SOMETHING WRONG, this->var_length[" <<v_index<<"]= "<<this->var_length[v_index]<<"  loop:  field_index= "<< field_index);
	      }
	      //}

	      if(i==0 && k==0)
	      {
		  vtkDebugMacro(<< "interpolateAndCopyContinuumData: my_rank= " << my_rank<<" : end v_index= " <<v_index<<" loop:  field_index= "<< field_index);
	      }

// 	      if(this->GetPointArrayStatus(v_index) == 1)
// 	      {
// 		  if(this->var_length[v_index] == 1)
// 		  {
// 		      scalars[cur_scalar_index++]->SetValue(index, data_array[field_index++][i]);
// 		  }
// 		  else if(this->var_length[v_index] == 3)
// 		  {
// 		      vectors[cur_vector_index++]->SetTuple3(index, 
// 		      					     data_array[field_index++][i],
// 		      					     data_array[field_index++][i],
// 		      					     data_array[field_index++][i]);
// 		  }
// 	      }


//         if(this->GetPointArrayStatus("Velocity"))
//           {
//           vectors->SetTuple3(
//             index,
//             data_array[0][i],
//             data_array[1][i],
//             data_array[2][i]);
//           }
//         if(this->GetPointArrayStatus("Pressure"))
//           {
//           pressure->SetValue(index, data_array[3][i]);
//           }
//         for(int kk=0; kk<this->num_vars; kk++)
//           {
//           // var_num is the index into the data_array, which is assumed to have
//           // uvwp (4) + extra vars.
//           int var_num = kk+4;
//           if(this->GetPointArrayStatus(this->var_names[kk]))
//             {
//             scalars[kk]->SetValue(index, data_array[var_num][i]);
//             }
//           }
	  } // for(v_index=0; v_index< this->num_vars; v_index++)
	  //}
        index ++;
	  
        }// for(i = 0; i < ntot; ++i)
      }
    }//for(k = 0; k < this->master[0]->nel; ++k)

  // We have now added all of the values to needed arrays,
  // now add those arrays to the UGrid if needed, 
  //   or remove if present and not needed

  vtkDebugMacro(<< "interpolateAndCopyContinuumData: my_rank= " << my_rank<<": we've added values to arrays, add them to ugrids, if needed");
  cur_scalar_index=0;
  cur_vector_index=0;
  // for each variable
  for(v_index=0; v_index< this->num_vars; v_index++)
  {
      // if it is requested
      if(this->GetPointArrayStatus(v_index) == 1)
      {
	  // if it is a scalar
	  if(this->var_length[v_index] == 1)
	  {
	      vtkDebugMacro(<< "interpolateAndCopyContinuumData: my_rank= " << my_rank<<": var["<<v_index<<"]: add scalars["<<cur_scalar_index<<"]");
	      this->UGrid->GetPointData()->AddArray(scalars[cur_scalar_index++]);
	  }
	  // if it is a vector 
	  else if(this->var_length[v_index] == 3)
	  {
	      vtkDebugMacro(<< "interpolateAndCopyContinuumData: my_rank= " << my_rank<<": var["<<v_index<<"]: add vectors["<<cur_vector_index<<"]");
	      this->UGrid->GetPointData()->AddArray(vectors[cur_vector_index++]);
	  }

      }// if(this->GetPointArrayStatus(v_index) == 1)
      else // User doesn't want it, remove it if already have it
      {
	  vtkDebugMacro(<< "interpolateAndCopyContinuumData: my_rank= " << my_rank<<": var["<<v_index<<"]: User doesn't want: "<<this->var_names[v_index]<<" : try to remove");
          // Does PV already have this array?  If so, remove it.
	  if (pv_ugrid->GetPointData()->GetArray(this->var_names[v_index]) != NULL)
	  {
	      pv_ugrid->GetPointData()->RemoveArray(this->var_names[v_index]);
	  }
	  // Do I already have this array?  If so, remove it.
	  if (this->UGrid->GetPointData()->GetArray(this->var_names[v_index]) != NULL)
	  {
	      this->UGrid->GetPointData()->RemoveArray(this->var_names[v_index]);
	  }
      }// else 
  }// for(v_index=0; v_index< this->num_vars; v_index++)

  // now delete the float arrays, since we're done with them

  for(cur_scalar_index=0; cur_scalar_index < this->num_used_scalars; cur_scalar_index++)
  {
      scalars[cur_scalar_index]->Delete();
  }
  if(scalars!=NULL)
      free(scalars);
  for(cur_vector_index=0; cur_vector_index < this->num_used_vectors; cur_vector_index++)
  {
      vectors[cur_vector_index]->Delete();
  }
  if(vectors!=NULL)
      free(vectors);



//   if(this->GetPointArrayStatus("Pressure"))
//     {
//     this->UGrid->GetPointData()->SetScalars(pressure);
//     }
//   else  // user does not want this variable
//     {
//     // Does PV already have this array?  If so, remove it.
//     if (pv_ugrid->GetPointData()->GetArray("Pressure") != NULL)
//       {
//       pv_ugrid->GetPointData()->RemoveArray("Pressure");
//       }
//     // Do I already have this array?  If so, remove it.
//     if (this->UGrid->GetPointData()->GetArray("Pressure") != NULL)
//       {
//       this->UGrid->GetPointData()->RemoveArray("Pressure");
//       }
//     }
//   pressure->Delete();

//   for(int kk=0; kk<this->num_vars; kk++)
//     {
//     if(this->GetPointArrayStatus(this->var_names[kk]))
//       {
//       this->UGrid->GetPointData()->AddArray(scalars[kk]);
//       }
//     else  // user does not want this variable
//       {
//       // Does PV already have this array?  If so, remove it.
//       if (pv_ugrid->GetPointData()->GetArray(this->var_names[kk]) != NULL)
//         {
//         pv_ugrid->GetPointData()->RemoveArray(this->var_names[kk]);
//         }
//       // Do I already have this array?  If so, remove it.
//       if (this->UGrid->GetPointData()->GetArray(this->var_names[kk]) != NULL)
//         {
//         this->UGrid->GetPointData()->RemoveArray(this->var_names[kk]);
//         }
//       }
//     scalars[kk]->Delete();
//     }

//   if(this->GetPointArrayStatus("Velocity"))
//     {
//     this->UGrid->GetPointData()->AddArray(vectors);
//     }
//   else  // user does not want this variable
//     {
//     // Does PV already have this array?  If so, remove it.
//     if (pv_ugrid->GetPointData()->GetArray("Velocity") != NULL)
//       {
//       pv_ugrid->GetPointData()->RemoveArray("Velocity");
//       }
//     // Do I already have this array?  If so, remove it.
//     if (this->UGrid->GetPointData()->GetArray("Velocity") != NULL)
//       {
//       this->UGrid->GetPointData()->RemoveArray("Velocity");
//       }
//     }

//   vectors->Delete();




  vtkDebugMacro(<< "interpolateAndCopyContinuumData: my_rank= " << my_rank<<"  EXIT");

}// vtkNektarReader::interpolateAndCopyContinuumData()

void vtkNektarReader::interpolateAndCopyBoundaryPoints(int alloc_res, int interp_res, vtkPoints* boundary_points)
{
  int my_rank = vtkMultiProcessController::GetGlobalController()->GetLocalProcessId();
  vtkDebugMacro(<< "interpolateAndCopyBoundaryPoints: my_rank= " << my_rank<<": ENTER: alloc_res= "<< alloc_res<< " interp_res: " << interp_res);
  int index=0;
  double *wk = dvector(0,alloc_res*alloc_res);
  Coord    boundary_X;

  vtkDebugMacro(<<"vtkNektarReader::interpolateAndCopyBoundaryPoints: alloc_res = " << alloc_res << ":  interp_res = " << interp_res);
  boundary_X.x = dvector(0,alloc_res*alloc_res*alloc_res-1);
  boundary_X.y = dvector(0,alloc_res*alloc_res*alloc_res-1);
  boundary_X.z = dvector(0,alloc_res*alloc_res*alloc_res-1);

  Tri tri;
  Quad quad;
  Nek_Facet_Type element_type = this->master[0]->flist[0]->identify();

  int ntot;
  vtkDebugMacro(<< "vtkNektarReader::vtkNektarReader():  tri->Nverts= "<< tri.Nverts);
  vtkDebugMacro(<< "vtkNektarReader::vtkNektarReader():  tri->dim()= "<< tri.dim());

  for(Bndry* B=this->Ubc; B; B = B->next)
    {
    // if it is a wall element
    if(B->type == 'W' || B->type == 'L' || B->type == 's' || B->type == 'G' || B->type == 'C' || B->type == 'A')
      {
      Element* F  = B->elmt;
#if 1
      switch(element_type)
      {
        case Nek_Tet: case Nek_Prism:
        {
        tri.qa = F->qa;
        tri.qb = F->qc; /* fix for prisms ?? */
        F->coord(&boundary_X);
        F->GetFace(boundary_X.x,B->face,wk);
        F->InterpToFace1(B->face,wk,boundary_X.x);
        Interp_symmpts(&tri,interp_res,boundary_X.x,boundary_X.x,'n'); //wk has coordinates of "x"
        F->GetFace(boundary_X.y,B->face,wk);
        F->InterpToFace1(B->face,wk,boundary_X.y);
        Interp_symmpts(&tri,interp_res,boundary_X.y,boundary_X.y,'n'); //wk has coordinates of "y"
        F->GetFace(boundary_X.z,B->face,wk);
        F->InterpToFace1(B->face,wk,boundary_X.z);
        ntot = Interp_symmpts(&tri,interp_res,boundary_X.z,boundary_X.z,'n'); //wk has coordinates of "z"
        break;
        }
        case Nek_Hex:
        {
        quad.qa = F->qa;
        quad.qb = F->qb;
        F->coord(&boundary_X);
        F->GetFace(boundary_X.x,B->face,wk);
        F->InterpToFace1(B->face,wk,boundary_X.x);
        Interp_symmpts(&quad,interp_res,boundary_X.x,boundary_X.x,'n'); //wk has coordinates of "x"
        F->GetFace(boundary_X.y,B->face,wk);
        F->InterpToFace1(B->face,wk,boundary_X.y);
        Interp_symmpts(&quad,interp_res,boundary_X.y,boundary_X.y,'n'); //wk has coordinates of "y"
        F->GetFace(boundary_X.z,B->face,wk);
        F->InterpToFace1(B->face,wk,boundary_X.z);
        ntot = Interp_symmpts(&quad,interp_res,boundary_X.z,boundary_X.z,'n'); //wk has coordinates of "z"
        break;
        }
      }
#else
      //vtkDebugMacro(<<  "interpolateAndCopyBoundaryPoints: my_rank= " << my_rank<<": Wall eleme: (B->face)==3: post assign tri.qb");
      // get the coordinates for this element
      F->coord(&boundary_X);

      F->GetFace(boundary_X.x,B->face,wk);
      F->InterpToFace1(B->face,wk,boundary_X.x);
      Interp_symmpts(&tri,interp_res,boundary_X.x,boundary_X.x,'n'); //wk has coordinates of "x"
      F->GetFace(boundary_X.y,B->face,wk);
      F->InterpToFace1(B->face,wk,boundary_X.y);
      Interp_symmpts(&tri,interp_res,boundary_X.y,boundary_X.y,'n'); //wk has coordinates of "y"
      F->GetFace(boundary_X.z,B->face,wk);
      F->InterpToFace1(B->face,wk,boundary_X.z);
      int ntot = Interp_symmpts(&tri,interp_res,boundary_X.z,boundary_X.z,'n'); //wk has coordinates of "z"
#endif

      // put all of the points for this boundary element into the list of boundary points
      for(int i = 0; i < ntot; ++i)
        {
        boundary_points->InsertPoint(index, boundary_X.x[i], boundary_X.y[i], boundary_X.z[i]);
        index ++;
        }
      }// if(B->type == 'W' || B->type == 'L' || B->type == 's' || B->type == 'G')
    }// for(B=this->Ubc; B; B = B->next)
  
  vtkDebugMacro(<< "interpolateAndCopyBoundaryPoints: my_rank= " << my_rank<<": num added= "<< index);


  free(boundary_X.x); 
  free(boundary_X.y); 
  free(boundary_X.z);
  free(wk);


  vtkDebugMacro(<< "interpolateAndCopyBoundaryPoints: my_rank= " << my_rank<<": EXIT ");

} // interpolateAndCopyBoundaryPoints()

void vtkNektarReader::interpolateAndCopyBoundaryData(int alloc_res, int num_verts, int interp_res)
{
  int my_rank = vtkMultiProcessController::GetGlobalController()->GetLocalProcessId();
  vtkDebugMacro(<< "interpolateAndCopyBoundaryData(): my_rank= " << my_rank<<"  ENTER");

  int ntot =3;
  int index=0;

  vtkDebugMacro(<< "interpolateAndCopyBoundaryData(): my_rank= " << my_rank<<" :: num_vert= "<< num_verts<<" : interp_res= " <<interp_res << "  : (interp_res+1)*interp_res/2= " << ((interp_res+1)*interp_res/2));

  double * boundary_valX = dvector(0, alloc_res*alloc_res*alloc_res-1);
  double * boundary_valY = dvector(0, alloc_res*alloc_res*alloc_res-1);
  double * boundary_valZ = dvector(0, alloc_res*alloc_res*alloc_res-1);

  int boundary_val_offset =0;

  vtkFloatArray* wall_shear_stress;
  Tri tri;
  Nek_Facet_Type element_type = this->master[0]->flist[0]->identify();
 
  /* Paris: Need to update this for quads*/
  if(this->GetDerivedVariableArrayStatus("Wall Shear Stress"))
    {
    wall_shear_stress = vtkFloatArray::New();
    wall_shear_stress->SetNumberOfComponents(3);
    wall_shear_stress->SetNumberOfTuples(num_verts);
    wall_shear_stress->SetName("Wall Shear Stress");

    // for each element in the boundary
    for(Bndry* B=this->Ubc; B; B = B->next)
      {
      Element* F  = B->elmt;
      // if it is a wall element
      if(B->type == 'W' || B->type == 'L' || B->type == 's' || B->type == 'G' || B->type == 'C' || B->type == 'A')
        {
        if(element_type == Nek_Tet)
        {
        tri.qa = F->qa;
        tri.qb = F->qc; /* fix for prisms */ /* Paris: need to check this*/
        }
        else if (element_type == Nek_Prism){
          if(B->face == 2 || B->face == 4){
            tri.qa = F->qa;
            tri.qb = F->qc;
          }
          else
          {
          tri.qa = F->qa;
          tri.qb = F->qb;
          }
        }
        else
        {
        tri.qa = F->qa;
        tri.qb = F->qb;
        }
#if 0
        Element* F  = B->elmt;
        //fprintf(stderr, "F->Nfverts(B->face): %d\n", F->Nfverts(B->face));

        if(F->Nfverts(B->face) == 3)
          {
          tri.qa = F->qa;
          tri.qb = F->qc; /* fix for prisms */
          }
        else
          {
          tri.qa = F->qa;
          tri.qb = F->qb;
          }
#endif

        ntot = Interp_symmpts(&tri,interp_res,this->WSS_all_vals[0]+boundary_val_offset, boundary_valX, 'n');
        ntot = Interp_symmpts(&tri,interp_res,this->WSS_all_vals[1]+boundary_val_offset, boundary_valY, 'n');
        ntot = Interp_symmpts(&tri,interp_res,this->WSS_all_vals[2]+boundary_val_offset, boundary_valZ, 'n');

        // put all of the points for this boundary element into the list of boundary points
        for(int i = 0; i < ntot; ++i)
          {
	      wall_shear_stress->SetTuple3(index, -boundary_valX[i], -boundary_valY[i], -boundary_valZ[i]);
	      index ++;
          }

        boundary_val_offset += tri.qa*tri.qb;
        }// if(B->type == 'W' || B->type == 'L' || B->type == 's' || B->type == 'G')

      }// for(B=this->Ubc; B; B = B->next)

    if(this->Boundary_UGrid)
      {
      this->Boundary_UGrid->GetPointData()->AddArray(wall_shear_stress);
      }
    else
      {
      vtkDebugMacro(<< "interpolateAndCopyBoundaryData(): my_rank= " << my_rank<<" boundary_ugrid is NULL");
      }

    wall_shear_stress->Delete();

    vtkDebugMacro(<< "~~~~~~~~~~~~~~~  interpolateAndCopyBoundaryData: my_rank= " << my_rank<<": num_verts= "<<num_verts<<"  index= "<<index);
    vtkDebugMacro(<< "~~~~~~~~~~~~~~~  interpolateAndCopyBoundaryData: my_rank= " << my_rank<<": ntot= "<<ntot<<"  (interp_res+1)*interp_res/2= "<<((interp_res+1)*interp_res/2));
    }// if(this->GetDerivedVariableArrayStatus("Wall Shear Stress"))

  free(boundary_valX);
  free(boundary_valY);
  free(boundary_valZ);
  vtkDebugMacro(<< "interpolateAndCopyBoundaryData(): my_rank= " << my_rank<<"  EXIT");

} // interpolateAndCopyBoundaryData()


// see if the current object is missing data that was requested
// return true if it is, otherwise false
bool vtkNektarReader::isObjectMissingData()
{
	int my_rank = vtkMultiProcessController::GetGlobalController()->GetLocalProcessId();
	
   if(this->curObj->resolution != this->ElementResolution)
   {
	   return(true);
   }
   vtkDebugMacro(<< "isObjectMissingData(): my_rank= " << my_rank<<" :  this->dynamic_coord_index = "<< this->dynamic_coord_index<< "  this->curObj->dynamic_mesh= "<<this->curObj->dynamic_mesh <<" this->DynamicMesh= "<<this->DynamicMesh);
   //if we can do dynamic mesh, AND: it has been turned on or off, OR it is on AND scale has changed
   if(this->dynamic_coord_index >=0 && ((this->curObj->dynamic_mesh != this->DynamicMesh) || (this->DynamicMesh && (this->curObj->dynamic_mesh_scale != this->DynamicMeshScale)) ) )
   {
	   vtkDebugMacro(<< "isObjectMissingData(): my_rank= " << my_rank<<" : mesh test failed, should return true");
	   return(true);	   
   }
   if(this->curObj->boundary_resolution != this->BoundaryResolution)
   {
	   return(true);
   }
   if(this->curObj->use_projection != this->UseProjection)
   {
	   return(true);	   
   }
   // check the stored variables
   for(int i=0; i<this->num_vars; i++)
   {
	   if(this->GetPointArrayStatus(i) ==1 && !this->curObj->vars[i])
	   {
		   return(true);
	   }
   }
   // check derived variables 
   if((this->GetDerivedVariableArrayStatus("Vorticity") && !this->curObj->vorticity) ||
      (this->GetDerivedVariableArrayStatus("lambda_2") && !this->curObj->lambda_2) ||
      (this->GetDerivedVariableArrayStatus("Wall Shear Stress") && !this->curObj->wss) ||
      (this->GetDerivedVariableArrayStatus("Stress Tensor") && !this->curObj->stress_tensor))
   {
	   return(true);
   }

   // object is not missing data
   return(false);
	
}// vtkNektarReader::isObjectMissingData()

bool vtkNektarReader::objectMatchesRequest()
{
// see if the current object matches the requested data
// return false if it does not match, otherwise true	
    int my_rank = vtkMultiProcessController::GetGlobalController()->GetLocalProcessId();

    if(curObj->resolution != this->ElementResolution)
    {
	    return(false);
    }
    if(curObj->boundary_resolution != this->BoundaryResolution)
    {
	    return(false);
    }
    if(this->dynamic_coord_index >=0 && ((this->curObj->dynamic_mesh != this->DynamicMesh) || (this->curObj->dynamic_mesh_scale != this->DynamicMeshScale )) )
    {
	    vtkDebugMacro(<< "objectMatchesRequest(): my_rank= " << my_rank<<" : mesh test failed, should return false");
	    return(false);	   
    }
    if(curObj->use_projection != this->UseProjection)
    {
	    return(false);
    }
    // check the stored variables
    for(int i=0; i<this->num_vars; i++)
    {
	    if(this->GetPointArrayStatus(i) != this->curObj->vars[i])
	    {
		    return(false);
	    }
    }
    // check derived variables
    if(this->GetDerivedVariableArrayStatus("Vorticity")  != this->curObj->vorticity)
    {
	    return(false);
    }
    if(this->GetDerivedVariableArrayStatus("lambda_2") != this->curObj->lambda_2)
    {
	    return(false);
    }	    
    if(this->GetDerivedVariableArrayStatus("Wall Shear Stress") != this->curObj->wss )
    {
	    return(false);
    }
    if(this->GetDerivedVariableArrayStatus("Stress Tensor") != this->curObj->stress_tensor )
    {
	    return(false);
    }

    // if we have the boundary geometry, but the current object has no bounday ugrid, or it has no cells
    
    if(this->HAVE_BOUNDARY_GEOM_FLAG)
    {
	    if(NULL == this->curObj->boundary_ugrid)
	    {
		    return(false);
	    }
	    if(this->curObj->boundary_ugrid->GetNumberOfCells()<=0)
	    {
		    return(false);
	    }
    }
    return(true);
	
}// vtkNektarReader::objectMatchesRequest()


bool vtkNektarReader::objectHasExtraData()
{
// see if the current object has extra data than was requested
// return false if object has less than request, otherwise true
	
        int my_rank = vtkMultiProcessController::GetGlobalController()->GetLocalProcessId();

	// if either resolution is different
	if((this->curObj->resolution != this->ElementResolution) ||
	   (this->curObj->boundary_resolution != this->BoundaryResolution))
	{
		return(false);
	}
	if(this->dynamic_coord_index >=0 && (this->curObj->dynamic_mesh != this->DynamicMesh))
	{
		return(false);	   
	}
	if(this->dynamic_coord_index >=0 && ((this->DynamicMesh) && (this->curObj->dynamic_mesh_scale != this->DynamicMeshScale )) )
	{
	    vtkDebugMacro(<< "objectHasExtraData(): my_rank= " << my_rank<<" : mesh test failed, should return false");
	    return(false);	   
	}
	// if projection does not match
	if(this->curObj->use_projection != this->UseProjection)
	{
		return(false);
	}
	// check the stored variables, if it was requested, but it is not in the current object, return false
	for(int i=0; i<this->num_vars; i++)
	{
		if(this->GetPointArrayStatus(i) && !this->curObj->vars[i])
		{
			return(false);
		}
	}
	
        // check derived variables, if it was requested, but it is not in the current object, return false
	if(this->GetDerivedVariableArrayStatus("Vorticity") && !this->curObj->vorticity)
	{
		return(false);
	}		
	if(this->GetDerivedVariableArrayStatus("lambda_2") && !this->curObj->lambda_2)
	{
		return(false);
	}		
	if(this->GetDerivedVariableArrayStatus("Wall Shear Stress") && !this->curObj->wss)
	{
		return(false);
	}
	if(this->GetDerivedVariableArrayStatus("Stress Tensor") && !this->curObj->stress_tensor)
	{
		return(false);
	}

	// if we have the boundary geometry, but the current object has no bounday ugrid, or it has no cells

	if(this->HAVE_BOUNDARY_GEOM_FLAG)
	{
		if(NULL == this->curObj->boundary_ugrid)
		{
			return(false);
		}
		if(this->curObj->boundary_ugrid->GetNumberOfCells()<=0)
		{
			return(false);
		}
	}

	vtkDebugMacro(<< "objectHasExtraData(): my_rank= " << my_rank<<" : returning true");
	return(true);
		
}// vtkNektarReader::objectHasExtraData()
