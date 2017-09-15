#ifndef STELLA_VTKWRITER_H
#define STELLA_VTKWRITER_H

#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTrivialProducer.h>
#include <vtkNew.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLPStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>

#include <stella/grid.h>

namespace stella { namespace io {

void writevtk(std::string fname, grid<2> & sgrid, cedar::cdr2::mpi::grid_func & gf)
{
	using namespace cedar;
	vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkDoubleArray> vals = vtkSmartPointer<vtkDoubleArray>::New();

	vals->SetNumberOfComponents(1);

	len_t nlx = sgrid.nlocal(0) - 1;
	len_t nly = sgrid.nlocal(1) - 1;
	len_t iend = sgrid.ibound(0) + 1;
	len_t jend = sgrid.ibound(1) + 1;
	len_t iext = sgrid.is(0) + sgrid.nlocal(0) - 3;
	len_t jext = sgrid.is(1) + sgrid.nlocal(1) - 3;

	if (sgrid.coord(0) == sgrid.nproc(0) - 1) {
		nlx--; iend--; iext--;
	}
	if (sgrid.coord(1) == sgrid.nproc(1) - 1) {
		nly--; jend--; jext--;
	}

	vals->SetNumberOfTuples(nlx * nly);
	vals->SetName("test");

	len_t acc = 0;
	for (len_t j = sgrid.ibase(1); j < jend; j++) {
		for (len_t i = sgrid.ibase(0); i < iend; i++) {
			points->InsertNextPoint(sgrid.vertex(0)(i,j), sgrid.vertex(1)(i,j), 1);
			vals->SetValue(acc, gf(i,j));
			acc++;
		}
	}

	structuredGrid->SetExtent(
		sgrid.is(0) - 1, iext,
		sgrid.is(1) - 1, jext,
		0, 0);

	structuredGrid->SetPoints(points);
	structuredGrid->GetPointData()->AddArray(vals);

	vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
		vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
	std::string fname_piece(fname + "_" + std::to_string(sgrid.mpi_rank()) + ".vts");
	writer->SetFileName(fname_piece.c_str());

	#if VTK_MAJOR_VERSION <= 5
	writer->SetInput(structuredGrid);
	#else
	writer->SetInputData(structuredGrid);
	#endif
	writer->Write();
}

}}


#endif
