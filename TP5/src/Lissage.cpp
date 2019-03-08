#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>
#include <fstream>
#include <queue>
#include <set>
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
	file << "COFF" << std::endl << mesh.size_of_vertices() << ' ' << mesh.size_of_facets() << " 0" << std::endl;
    for (Vertex_iterator i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i) {
		file << coords[i];
		file<<' '<< (int)((sqrt(CGAL::squared_distance(coords[i], i->point())) - delta_min)/(delta_max-delta_min) * 255 )<<' '<<"250"<<' '<<"250" <<' '<<"1.0"<< std::endl;
		//file<<' '<< "200" <<' '<<"200"<<' '<<"200" <<' '<<"1.0"<< std::endl;
	}
    for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i) {
        Halfedge_facet_circulator j = i->facet_begin();
        CGAL_assertion( CGAL::circulator_size(j) >= 3);
        file << CGAL::circulator_size(j) << ' ';
        do {
            file << ' ' << std::distance(mesh.vertices_begin(), j->vertex());
        } while ( ++j != i->facet_begin());
		file << std::endl;
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
	float sum_x = 0;
	float sum_y = 0;
	float sum_z = 0;
	int nbVoisins;
	float delta;
	delta_min = FLT_MAX;
	delta_max = 0;
	
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
		delta = sqrt(CGAL::squared_distance(new_p, coords[i]));
		if(delta>delta_max)
			delta_max = delta;
		if(delta<delta_min)
			delta_min = delta;
		sommets[i] = new_p;
	}
	return sommets;
}

Vertex_coord_map laplacianSmoothingRadius(Polyhedron & mesh, Vertex_coord_map & coords, double radius){
	Vertex_coord_map sommets;
	float sum_x, sum_y, sum_z;
	int nbVoisins;
	float delta;
	delta_min = FLT_MAX;
	delta_max = -1;
	std::set<Vertex_handle> visite;
	std::queue<Vertex_handle> aVisiter;

	for (Vertex_iterator i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i) {
		nbVoisins = 0;
		sum_x = 0; sum_y = 0; sum_z = 0;
		aVisiter.push(i);
		visite.insert(i);
		while(!aVisiter.empty()){
			Vertex_handle vertCour = aVisiter.front();
			aVisiter.pop();
			sum_x += coords[vertCour].x();
			sum_y += coords[vertCour].y();
			sum_z += coords[vertCour].z();
			nbVoisins ++;
			Polyhedron::Halfedge_around_vertex_circulator e(vertCour->vertex_begin()), e_end(e);
			do{
				Vertex_handle new_v = e->opposite()->vertex();
				if(visite.find(new_v)==visite.end() && CGAL::squared_distance(coords[new_v], coords[i]) < radius*radius){
					aVisiter.push(new_v);
					visite.insert(new_v);
				}
			} while (++e != e_end);
		}
		visite.clear();
		Point_3 new_p = Point_3 (sum_x/(float)nbVoisins, sum_y/(float)nbVoisins, sum_z/(float)nbVoisins);
		delta = sqrt(CGAL::squared_distance(new_p, coords[i]));
		if(delta>delta_max)
			delta_max = delta;
		if(delta<delta_min)
			delta_min = delta;
		sommets[i] = new_p;
	}
	return sommets;
}

Vertex_coord_map laplacianSmoothingRadiusPond(Polyhedron & mesh, Vertex_coord_map & coords, double radius){
	Vertex_coord_map sommets;
	float sum_x, sum_y, sum_z;
	int nbVoisins;
	float delta;
	delta_min = FLT_MAX;
	delta_max = -1;
	std::set<Vertex_handle> visite;
	std::queue<Vertex_handle> aVisiter;

	for (Vertex_iterator i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i) {
		nbVoisins = 0;
		sum_x = 0; sum_y = 0; sum_z = 0;
		aVisiter.push(i);
		visite.insert(i);
		while(!aVisiter.empty()){
			Vertex_handle vertCour = aVisiter.front();
			aVisiter.pop();
			int poids = (int)((1 - (sqrt(squared_distance(coords[vertCour], coords[i]))/radius))*10);
			sum_x += coords[vertCour].x() * poids;
			sum_y += coords[vertCour].y() * poids;
			sum_z += coords[vertCour].z() * poids;
			nbVoisins += poids;
			Polyhedron::Halfedge_around_vertex_circulator e(vertCour->vertex_begin()), e_end(e);
			do{
				Vertex_handle new_v = e->opposite()->vertex();
				if(visite.find(new_v)==visite.end() && CGAL::squared_distance(coords[new_v], coords[i]) < radius*radius){
					aVisiter.push(new_v);
					visite.insert(new_v);
				}
			} while (++e != e_end);
		}
		visite.clear();
		Point_3 new_p = Point_3 (sum_x/(float)nbVoisins, sum_y/(float)nbVoisins, sum_z/(float)nbVoisins);
		delta = sqrt(CGAL::squared_distance(new_p, coords[i]));
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

	Vertex_coord_map sommetsLaplacien;
	//sommetsLaplacien = laplacianSmoothing(mesh, sommetsInitiaux);
	//sommetsLaplacien = laplacianSmoothingRadius(mesh, sommetsInitiaux, 1.5);
	sommetsLaplacien = laplacianSmoothingRadiusPond(mesh, sommetsInitiaux, 1.5);

	std::cout << "delta min: " << delta_min << std::endl;
	std::cout << "delta max: " << delta_max << std::endl;
	draw(mesh,sommetsLaplacien);
	
	return 0;
}
