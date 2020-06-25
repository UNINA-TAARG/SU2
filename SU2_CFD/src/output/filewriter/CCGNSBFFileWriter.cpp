/*!
 * \file CCGNSBFFileWriter.cpp
 * \brief Filewriter class for Tecplot ASCII format.
 * \author T. Albring
 * \version 7.0.5 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../../include/output/filewriter/CCGNSBFFileWriter.hpp"
#include "../../../../externals/cgns/cgnslib.h"
#include "../../../../Common/include/geometry/CGeometry.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../../include/solvers/CSolver.hpp"
#include "../../../../Common/include/option_structure.hpp"
#include "../../../include/output/COutputLegacy.hpp"

#include <algorithm>

const string CCGNSBFFileWriter::fileExt = ".cgns";

CCGNSBFFileWriter::CCGNSBFFileWriter(string valFileName, CParallelDataSorter *valDataSorter) :
  CFileWriter(std::move(valFileName), valDataSorter, fileExt){}

CCGNSBFFileWriter::~CCGNSBFFileWriter(){}

void CCGNSBFFileWriter::Write_Data_BF(CConfig *config, CGeometry *geometry){

  if (!dataSorter->GetConnectivitySorted()){
    SU2_MPI::Error("Connectivity must be sorted.", CURRENT_FUNCTION);
  }

 #ifdef HAVE_CGNS

	/*	Internal Variables	*/

	int cgns_err, cgns_file, cgns_file1, cgns_flow, cgns_field, nbndry, NormalIndex, bc_wall, bc_far, bc_nulli, coord_index;
	int nbases, base_number, cell_dim, phys_dim, nzones, zone_number, ncoords, nsections, i, S;
	int iMarker, nbocos, ndataset, bcnum, parent_flag, rank;
	cgsize_t isize[3][1];
	cgsize_t start, end, ElementDataSize, NormalListSize, npnts, j, k, z, p, q, NewElementDataSize, dimensions;
	cgsize_t pnts[1][2];
	cgsize_t A[4], B[4];
	cgsize_t *Elements = NULL;
	cgsize_t *NewElements = NULL;
	PointSetType_t ptset_type;
	BCType_t bocotype;
	DataType_t datatype;
	GridLocation_t location;
	ElementType_t type;
	DataClass_t dataclass;
	su2double *rho = nullptr;
	su2double *rhou = nullptr;
	su2double *rhov = nullptr;
	su2double *rhow = nullptr;
	su2double *u = nullptr;
	su2double *v = nullptr;
	su2double *w = nullptr;
	su2double RefDensity, RefPressure, RefTemperature, RefVelocity;
    unsigned short iDim = 0, nDim = dataSorter->GetnDim();
	char boconame;
	char basename[] = "Base";
	char zonename[] = "Zone 1";
	char coordname[32];
	char ElementSectionName[32];
	unsigned long jVar;
	unsigned short NVar, iVar, iPoint;
	bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
	bool nondim1 = (config->GetRef_NonDim() == FREESTREAM_PRESS_EQ_ONE);
	bool nondim2 = (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_MACH);
	bool nondim3 = (config->GetRef_NonDim() == FREESTREAM_VEL_EQ_ONE);
	stringstream name, secname, ss;
	string mesh_file, base_file, wall_string, far_string, null_string, str;
    const vector<string> fieldNames = dataSorter->GetFieldNames();
    unsigned long nPoint_Global = dataSorter->GetnPointsGlobal();
    nElem = geometry->GetnElem();
    unsigned short nVar = fieldNames.size();

    su2double Data[nVar][nPoint_Global] = {0.0};

	if (nondim1) {
		RefDensity = 1.0;
		RefPressure = 1.0;
		RefTemperature = 1.0;
		cout << "Non-Dimensional simulation (P=1.0, Rho=1.0, T=1.0 at the farfield)." << endl;
	}
	else if (nondim2) {
		RefVelocity = config->GetMach();
		RefDensity = 1.0;
		RefTemperature = 1.0;
		cout << "Non-Dimensional simulation (V=Mach, Rho=1.0, T=1.0 at the farfield)." << endl;
	}
	else if (nondim3) {
		RefVelocity = 1.0;
		RefDensity = 1.0;
		RefTemperature = 1.0;
		cout << "Non-Dimensional simulation (V=1.0, Rho=1.0, T=1.0 at the farfield)." << endl;
	}

    for (iPoint = 0; iPoint < dataSorter->GetnPoints(); iPoint++) {
            for (iVar = 0; iVar < nVar; iVar++)
            	Data[iVar][iPoint] = dataSorter->GetData(iVar, iPoint);
          }

	isize[0][0] = (cgsize_t)nPoint_Global;        // vertex size
	isize[1][0] = (cgsize_t)nElem;        // cell size
	isize[2][0] = 0;		//(cgsize_t)nSurf_Elem;            // boundary vertex size (zero if elements not sorted)

	rho = new su2double[nPoint_Global];
	rhou = new su2double[nPoint_Global];
	u = new su2double[nPoint_Global];
	for (i = 0; i < nPoint_Global; i++)	{
		rho[i] = Data[3][i];
		rhou[i] = Data[4][i];
		u[i] = rhou[i]/rho[i];
	}
	rhov = new su2double[nPoint_Global];
	v = new su2double[nPoint_Global];
	for (i = 0; i < nPoint_Global; i++)	{
		rhov[i] = Data[5][i];
		v[i] = rhov[i]/rho[i];
	}
	rhow = new su2double[nPoint_Global];
	w = new su2double[nPoint_Global];
	for (i = 0; i < nPoint_Global; i++)	{
		rhow[i] = Data[6][i];
		w[i] = rhow[i]/rho[i];
	}
	/*	Open Mesh File	*/
	mesh_file = config->GetMesh_FileName();
	cgns_err = cg_open(mesh_file.c_str(), CG_MODE_READ, &cgns_file);
	if(cgns_err) cg_error_print();
	nbases = 1;
	nzones = 1;

	/*      Create Solution File    */
//	base_file = config->GetSolution_FileName();
	base_file = fileName;
//	base_file = base_file.append(".cgns");
	cgns_err = cg_open(base_file.c_str(), CG_MODE_WRITE, &cgns_file1);
	if(cgns_err) cg_error_print;
	cout << "Soulution file opened" << endl;
	cell_dim = 3;
	phys_dim = 3;
	cgns_err = cg_base_write(cgns_file1, basename, cell_dim, phys_dim,  &base_number);
	if(cgns_err) cg_error_print;
	cgns_err = cg_zone_write(cgns_file1, base_number, zonename, *isize, Unstructured, &zone_number);

	/*	Write Coordinates	*/
	cgns_err = cg_ncoords(cgns_file, nbases, nzones, &ncoords);
	if(cgns_err) cg_error_print();
	for ( i = 1; i < ncoords + 1; i++)	{
		cgns_err = cg_coord_info(cgns_file, nbases, nzones, i, &datatype, coordname);
		if(cgns_err) cg_error_print();
		start = 1;
		end = isize[0][0];
		vector<double> buf(end);
		cgns_err = cg_coord_read(cgns_file, nbases, nzones, coordname,
				datatype, &start, &end, buf.data());
		if(cgns_err) cg_error_print();
		cgns_err = cg_coord_write(cgns_file1, base_number, zone_number, datatype, coordname, buf.data(), &coord_index);
		if(cgns_err) cg_error_print();
	}

	/*	Write Connectivity	*/
	cgns_err = cg_nsections(cgns_file, nbases, nzones, &nsections);
	if(cgns_err) cg_error_print();
	for ( i = 1; i < nsections + 1; i++)	{
		cgns_err = cg_section_read(cgns_file, nbases, nzones, i,
				ElementSectionName, &type, &start,
				&end, &nbndry, &parent_flag);
		if(cgns_err) cg_error_print();
		cgns_err = cg_ElementDataSize(cgns_file, nbases, nzones, i, &ElementDataSize);
		if(cgns_err) cg_error_print();
		Elements = new cgsize_t[ElementDataSize];
		cgns_err = cg_elements_read(cgns_file, nbases, nzones, i, Elements, NULL);
		if(cgns_err) cg_error_print();
		cgns_err = cg_section_write(cgns_file1, base_number, zone_number,
				ElementSectionName, type, start,
				end, nbndry, Elements, &S);
		if(cgns_err) cg_error_print();

		/*	Add Boundary Condition	*/
		for (j = 0; j <= config->GetnMarker_All(); j++) {
			if (ElementSectionName == config->GetMarker_All_TagBound(j)) {
				if (config->GetMarker_All_KindBC(j) == 26)      {
					bocotype = (BCType_t)20;
					ptset_type = (PointSetType_t)PointList;
					cgns_err = cg_boco_write(cgns_file1, base_number, zone_number, ElementSectionName, bocotype, ptset_type, ElementDataSize, Elements, &bcnum);
					if(cgns_err) cg_error_print();
					location = Vertex;
					cgns_err = cg_boco_gridlocation_write(cgns_file1, base_number, zone_number, bcnum, location);
					if(cgns_err) cg_error_print();
				}
				else if (config->GetMarker_All_KindBC(j) == 2)  {
					bocotype = (BCType_t)7;
					ptset_type = (PointSetType_t)PointList;
					cgns_err = cg_boco_write(cgns_file1, base_number, zone_number, ElementSectionName, bocotype, ptset_type, ElementDataSize, Elements, &bcnum);
					if(cgns_err) cg_error_print();
					location = Vertex;
					cgns_err = cg_boco_gridlocation_write(cgns_file1, base_number, zone_number, bcnum, location);
					if(cgns_err) cg_error_print();
				}
				else    {
					bocotype = (BCType_t)0;
					ptset_type = (PointSetType_t)PointList;
					cgns_err = cg_boco_write(cgns_file1, base_number, zone_number, ElementSectionName, bocotype, ptset_type, ElementDataSize, Elements, &bcnum);
					if(cgns_err) cg_error_print();
					location = Vertex;
					cgns_err = cg_boco_gridlocation_write(cgns_file1, base_number, zone_number, bcnum, location);
					if(cgns_err) cg_error_print();
				}
			}
		}
	}

	cout << "Grid rewritten in the output file." << endl;
	cout << "Added Boundary Conditions." << endl;

	/*      Close CGNS Mesh File    */
	cgns_err = cg_close(cgns_file);
	if(cgns_err) cg_error_print();

	/*	Write Solution	*/
	if(compressible) {
		cgns_err = cg_sol_write(cgns_file1, base_number, zone_number, (char *)"Solution", Vertex, &cgns_flow);
		if(cgns_err) cg_error_print();
		cgns_err = cg_goto(cgns_file1, nbases, "Zone_t", 1, "end");
		dataclass = NormalizedByUnknownDimensional;
		cgns_err = cg_dataclass_write(dataclass);
		nbases = base_number;
		nzones = zone_number;
		if(cgns_err) cg_error_print();
		switch (config->GetKind_Solver()) {
		case EULER:
			if (cell_dim == 2)      {
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"zero", Data[0], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"x", Data[1], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"y", Data[2], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Density", Data[3], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumX", Data[4], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumY", Data[5], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Energy", Data[6], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Pressure", Data[8], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Temperature", Data[9], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Mach", Data[10], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cp", Data[11], &cgns_field);
				if(cgns_err) cg_error_print();
			}
			else {
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"zero", Data[0], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"x", Data[1], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"y", Data[2], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"z", Data[3], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Density", Data[4], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumX", Data[5], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumY", Data[6], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumZ", Data[7], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Energy", Data[8], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Pressure", Data[9], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Temperature", Data[10], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Mach", Data[11], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cp", Data[12], &cgns_field);
				if(cgns_err) cg_error_print();
			}
			break;
		case NAVIER_STOKES:
			if(cell_dim == 2){
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"zero", Data[0], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"x", Data[1], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"y", Data[2], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Density", Data[3], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumX", Data[4], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumY", Data[5], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Energy", Data[6], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Pressure", Data[7], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Temperature", Data[8], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Mach", Data[9], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cp", Data[10], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Viscosity", Data[11], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cf_x", Data[12], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cf_y", Data[13], &cgns_field);
				if(cgns_err) cg_error_print();
			}
			else {
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"zero", Data[0], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"x", Data[1], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"y", Data[2], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"z", Data[3], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Density", Data[4], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumX", Data[5], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumY", Data[6], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumZ", Data[7], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Energy", Data[8], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Pressure", Data[9], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Temperature", Data[10], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Mach", Data[11], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cp", Data[12], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Viscosity", Data[13], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cf_x", Data[14], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cf_y", Data[15], &cgns_field);
				if(cgns_err) cg_error_print();
				cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cf_z", Data[16], &cgns_field);
				if(cgns_err) cg_error_print();
			}
			break;
		case RANS:
			switch (config->GetKind_Turb_Model()) {
			case SA:
				if(cell_dim == 2){
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"x", Data[0], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"y", Data[1], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Density", Data[2], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumX", Data[3], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumY", Data[4], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Energy", Data[5], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Nu_Tilde", Data[6], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Pressure", Data[7], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Temperature", Data[8], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Mach", Data[9], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cp", Data[10], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Viscosity", Data[11], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cf_x", Data[12], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cf_y", Data[13], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Enthalpy", Data[14], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"y+", Data[15], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"EddyViscosity", Data[16], &cgns_field);
					if(cgns_err) cg_error_print();
				}
				else {
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Density", rho, &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumX", rhou, &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumY", rhov, &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumZ", rhow, &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"EnergyStagnationDensity", Data[7], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Pressure", Data[9], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"ViscosityEddy", Data[19], &cgns_field);
					if(cgns_err) cg_error_print();
				}
				break;
			case SST:
				if(cell_dim ==2) {
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"x", Data[0], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"y", Data[1], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Density", Data[2], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumX", Data[3], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumY", Data[4], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Energy", Data[5], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"TurbulentKineticEnergy", Data[6], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"TurbulentDissipationRate", Data[7], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Nu_Tilde", Data[8], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Pressure", Data[9], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Temperature", Data[10], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Mach", Data[11], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cp", Data[12], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Viscosity", Data[13], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cf_x", Data[14], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cf_y", Data[15], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Enthalpy", Data[16], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"y+", Data[17], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"EddyViscosity", Data[18], &cgns_field);
					if(cgns_err) cg_error_print();
				}
				else {
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"x", Data[0], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"y", Data[1], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"z", Data[2], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Density", Data[3], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumX", Data[4], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumY", Data[5], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"MomentumZ", Data[6], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Energy", Data[7], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"TurbulentKineticEnergy", Data[8], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"TurbulentDissipationRate", Data[9], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Nu_Tilde", Data[10], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Pressure", Data[11], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Temperature", Data[12], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Mach", Data[13], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cp", Data[14], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Viscosity", Data[15], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cf_x", Data[16], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cf_y", Data[17], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Cf_z", Data[18], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"Enthalpy", Data[19], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"y+", Data[20], &cgns_field);
					if(cgns_err) cg_error_print();
					cgns_err = cg_field_write(cgns_file1, nbases, nzones, cgns_flow, RealDouble,"EddyViscosity", Data[21], &cgns_field);
					if(cgns_err) cg_error_print();
				}
				break;
			}
			break;
		}
	}
	cout << "Added the Soution Data." << endl;


	/*	Close CGNS Solution File	*/
	cgns_err = cg_close(cgns_file1);
	if(cgns_err) cg_error_print();

	cout << "CGNS input file for BreakForce succesfully created!" << endl;

#else // Not built with CGNS support

	cout << "CGNS file requested but SU2 was built without CGNS support. No file written" << "\n";

#endif

}


