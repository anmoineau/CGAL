#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>
#include <fstream>
#include <queue>
#include <CGAL/squared_distance_3.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Kernel::Point_3 Point_3;
typedef Polyhedron::Vertex_handle Vertex_handle;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;
typedef std::map<Polyhedron::Vertex_handle, Point_3> Vertex_coord_map;

float delta_min;
float delta_max;

void draw(Polyhedron & mesh,Vertex_coord_map & coords){
	std::fstream file;
	file.open("../Rendu.off");
	file << "OFF" << std::endl << mesh.size_of_vertices() << ' ' << mesh.size_of_facets() << " 0" << std::endl;
    for (Vertex_iterator i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i) {
		file << coords[i] << std::endl;
	}
    for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i) {
        Halfedge_facet_circulator j = i->facet_begin();
        CGAL_assertion( CGAL::circulator_size(j) >= 3);
        file << CGAL::circulator_size(j) << ' ';
        do {
            file << ' ' << std::distance(mesh.vertices_begin(), j->vertex());
        } while ( ++j != i->facet_begin());
		file<<"  "<<"1.000"<<' '<<"1.000"<<' '<<"1.000"<<' '<<"0.75"<< std::endl;
	}
}

Vertex_coord_map getCoords(Polyhedron & mesh){
	Vertex_coord_map sommets;
	for (Vertex_iterator i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i) {
		sommets[i] = i->point();
	}
	return sommets;
}

Vertex_coord_map laplacianSmoothing(Polyhedron & mesh, Vertex_coord_map & coords){
	Vertex_coord_map sommets;
	Vertex_coord_map::iterator it;
	float sum_x = 0;
	float sum_y = 0;
	float sum_z = 0;
	float delta;
	delta_min = 1000000;
	delta_max = 0;

	int nbVoisins;
	for (Vertex_iterator i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
	{
		nbVoisins = 0;
		sum_x = 0; sum_y = 0; sum_z = 0;
		Polyhedron::Halfedge_around_vertex_circulator e(i->vertex_begin()), e_end(e);
		do{
			Vertex_handle new_v = e->opposite()->vertex();
			sum_x += coords[new_v].x();
			sum_y += coords[new_v].y();
			sum_z += coords[new_v].z();
			nbVoisins ++;
		} while (++e != e_end);
		Point_3 new_p = Point_3 (sum_x/nbVoisins, sum_y/nbVoisins, sum_z/nbVoisins);
		delta = CGAL::squared_distance(new_p, coords[i]);
		if(delta>delta_max)
			delta_max = delta;
		if(delta<delta_min)
			delta_min = delta;
		sommets[i] = new_p;
	}
	return sommets;
}

int main(int argc, char* argv[])
{
	if (argc < 2) {
		std::cerr << "Il manque un paramètre au programme. Veuillez lui donner en entrée un nom de fichier au format off." << std::endl;
		return 1;
	}

	Polyhedron mesh;
	std::ifstream input(argv[1]);
	if (!input || !(input >> mesh) || mesh.is_empty()) {
		std::cerr << "Le fichier donné n'est pas un fichier off valide." << std::endl;
		return 1;
	}
	Vertex_coord_map sommetsInitiaux = getCoords(mesh);
	std::cout << "nombre de sommets: " << sommetsInitiaux.size() << std::endl; 
	Vertex_coord_map sommetsLaplacien = laplacianSmoothing(mesh, sommetsInitiaux);
	draw(mesh,sommetsLaplacien);
	std::cout << "delta min: " << delta_min << std::endl;
	std::cout << "delta max: " << delta_max << std::endl;
	return 0;
}
