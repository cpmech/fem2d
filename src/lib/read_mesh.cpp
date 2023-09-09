#include <clocale>
#include <cwchar>
#include <memory>
#include <vector>

#include "read_mesh.h"

const size_t MAX_LINE_WIDTH = 500;
const size_t MAX_KIND_WIDTH = 24;

size_t find_first_non_empty_character(size_t nmax, wchar_t line[]) {
    for (size_t i = 0; i < nmax; i++) {
        if (line[i] != L' ') {
            return i;
        }
    }
    return nmax;
}

std::unique_ptr<CoordinatesAndConnectivity> read_mesh(const std::string &filename) {

    setlocale(LC_ALL, "en_US.UTF-8"); // important to read UTF8 files

    FILE *f = fopen(filename.c_str(), "r");
    if (f == NULL) {
        throw "read_mesh: cannot open file";
    }

    wchar_t line[MAX_LINE_WIDTH];

    bool reading_header = true;
    bool reading_coordinates = false;
    bool reading_connectivity = false;

    int nread;
    size_t ndim, npoint, ncell;
    size_t id, att, a, b, c;
    double x, y;
    wchar_t kind[MAX_KIND_WIDTH];

    size_t counter_points = 0;
    size_t counter_cells = 0;

    std::vector<double> coordinates;
    std::vector<size_t> connectivity;
    size_t element_nnode;

    while (fgetws(line, MAX_LINE_WIDTH, f) != NULL) {
        if (line[0] == L'#' || line[0] == L'\n') {
            continue;
        }
        // if (find_first_non_empty_character(line_max,line)==line_max)

        if (reading_header) {

            nread = swscanf(line, L"%zu %zu %zu", &ndim, &npoint, &ncell);
            if (nread != 3) {
                fclose(f);
                throw "read_mesh cannot parse the dimensions: ndim npoint ncell";
            }
            if (ndim != 2) {
                fclose(f);
                throw "read_mesh works with ndim=2 only at this time";
            }
            // coo = CooMatrix::make_new(layout, m, max);
            reading_header = false;
            reading_coordinates = true;
            coordinates.resize(npoint * ndim);

        } else if (reading_coordinates) {

            nread = swscanf(line, L"%zu %lg %lg", &id, &x, &y);
            if (nread != 3) {
                fclose(f);
                throw "read_mesh cannot parse the coordinate: id x y";
            }
            if (id != counter_points) {
                fclose(f);
                throw "read_mesh requires that the id and index of points must equal each other";
            }
            coordinates[id * ndim] = x;
            coordinates[id * ndim + 1] = y;
            counter_points++;
            if (counter_points == npoint) {
                reading_coordinates = false;
                reading_connectivity = true;
            }

        } else if (reading_connectivity) {

            nread = swscanf(line, L"%zu %zu %24ls %zu %zu %zu", &id, &att, kind, &a, &b, &c);
            if (!(nread == 5 || nread == 6)) {
                fclose(f);
                throw "read_mesh cannot parse the connectivity: id att kind a b [c]";
            }
            if (id != counter_cells) {
                fclose(f);
                throw "read_mesh requires that the id and index of cells must equal each other";
            }
            if (wcsncmp(kind, L"lin2", 4) == 0) {
                element_nnode = 2;
            } else if (wcsncmp(kind, L"tri3", 4) == 0) {
                element_nnode = 3;
            } else {
                fclose(f);
                throw "read_mesh works with lin2 and lin3 only at this time";
            }
            if (element_nnode == 2 && nread != 5) {
                fclose(f);
                throw "read_mesh cannot read the lin2 cell connectivity";
            }
            if (element_nnode == 3 && nread != 6) {
                fclose(f);
                throw "read_mesh cannot read the tri3 cell connectivity";
            }
            if (connectivity.size() == 0) {
                connectivity.resize(ncell * element_nnode);
            }
            connectivity[id * element_nnode] = a;
            connectivity[id * element_nnode + 1] = b;
            if (element_nnode == 3) {
                connectivity[id * element_nnode + 2] = c;
            }
            counter_cells++;
            if (counter_cells == ncell) {
                break;
            }
        }
    }

    fclose(f);

    if (reading_header) {
        throw "read_mesh failed to read the header";
    }
    if (reading_coordinates) {
        throw "read_mesh failed to read the coordinates";
    }
    if (!reading_connectivity) {
        throw "read_mesh failed to read the connectivity";
    }
    if (counter_cells != ncell) {
        throw "read_mesh failed because there are not enough cell data";
    }

    return std::unique_ptr<CoordinatesAndConnectivity>{
        new CoordinatesAndConnectivity{
            coordinates,
            connectivity,
        }};
}