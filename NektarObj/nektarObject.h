#ifndef __nektarObject_h

#include "vtkUnstructuredGrid.h"

#define MAX_VARS 16

class nektarObject 
{
 public:
    vtkUnstructuredGrid* ugrid;
    vtkUnstructuredGrid* boundary_ugrid;
    bool vorticity;
    bool lambda_2;
    bool wss;
    bool stress_tensor;
    bool vars[MAX_VARS];
    bool der_vars[MAX_VARS];
    int index;
    int resolution;
    int boundary_resolution;
    bool use_projection;
    bool dynamic_mesh;
    double dynamic_mesh_scale;
    nektarObject *prev;
    nektarObject *next;
    char * dataFilename;

    nektarObject *find_obj(int id);
    void setDataFilename(char* filename);
    void reset();
   
// protected:
    nektarObject();
    ~nektarObject();

};


class nektarList
{
 public:
    nektarObject* head;
    nektarObject* tail;
    int max_count;
    int cur_count;
    nektarObject* getObject(int);

// protected:
    nektarList();
    ~nektarList();
};


#endif
