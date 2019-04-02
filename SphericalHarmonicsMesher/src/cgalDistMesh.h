#ifndef CGALDISTMESH_H
#define CGALDISTMESH_H

typedef struct {
        unsigned minNumPanels, maxNumPanels, verbose;
            double maxEdgeLength; // *initial* value, before acct'ing for minNumPanels, maxNumPanels
                bool convertToGMSH;
                    const char *offFile, *gmshFile;
} mesh_inputs;

typedef struct {
        unsigned numPanels;
            double surfArea, volume;
                double *panelAreas, *centroids;
                    const char *offFile, *gmshFile;
} mesh_outputs;

typedef double (*dist_fx)(double x, double y, double z);

void cgalDistMesh(dist_fx dist_function, double bound, mesh_inputs *mi, mesh_outputs *mo, double d_bound=INFINITY);
//void ConvertOFFtoGMSH(const char *OFFFileIn, const char *GMSHFileOut, int verbose=0);

#endif
