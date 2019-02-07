#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>
#include <fstream>
#include <CGAL/squared_distance_3.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Kernel::Point_3 Point_3;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;

void Draw(Polyhedron & P, std::map<Polyhedron::Facet_handle, int> categories){
	//liste de couleurs aleatoire de taille nombre de categorie -> Couleurs[categories[i]]
	std::fstream file;
	file.open("../Rendu.off");
	file << "OFF" << std::endl << P.size_of_vertices() << ' ' << P.size_of_facets() << " 0" << std::endl;
    std::copy( P.points_begin(), P.points_end(), std::ostream_iterator<Point_3>( file, "\n"));
    for (Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
        Halfedge_facet_circulator j = i->facet_begin();
        CGAL_assertion( CGAL::circulator_size(j) >= 3);
        file << CGAL::circulator_size(j) << ' ';
        do {
            file << ' ' << std::distance(P.vertices_begin(), j->vertex());
        } while ( ++j != i->facet_begin());
		file<<"  "<<((categories[i] == 1)? 0.25 : 0.75)<<' '<<((categories[i] == 2)? 1.000 : 0.000)<<' '<<"1.000"<<' '<<"0.75"<< std::endl;
        file << std::endl;
	}
}

std::map<Polyhedron::Facet_handle, double> ComputePerimetre(Polyhedron & P,float & min,float & max, float & mean){
	max = 0;
	min = 1000;
	float sum = 0;
	std::map<Polyhedron::Facet_handle, double> perimetres;
	for (Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		float perimetre = 0;
		Halfedge_facet_circulator j = i->facet_begin();
		CGAL_assertion( CGAL::circulator_size(j) >= 3);
		do {
            perimetre += CGAL::squared_distance(j->vertex()->point(),j->opposite()->vertex()->point());
        } while ( ++j != i->facet_begin());
		if(perimetre>max)
			max = perimetre;
		if(perimetre<min)
			min = perimetre;
		perimetres[i] = perimetre; 
		sum += perimetre;
	}
	mean = sum/perimetres.size();
	std::cout << "Perimetre moyen: " << mean << std::endl;
	return perimetres;
}

std::map<Polyhedron::Facet_handle, int> ComputeTreshold(Polyhedron & P, std::map<Polyhedron::Facet_handle, double> & perimetres, const float & mean){
	std::map<Polyhedron::Facet_handle, int> categories;

    for (Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		categories[i] = (perimetres[i]> mean)? 1 : 2;
	}
	return categories;
}

/*std::map<Polyhedron::Facet_handle, int> ModifyCategorie(Polyhedron & P, std::map<Polyhedron::Facet_handle, int> & categories){
	std::map<Polyhedron::Facet_handle, bool> visite;
	//parcours face -> visité
	//faire le tour des arretes -> empiler facet associée si non visitée 
	//attribuer nouvelle categorie
	//jusqu'a pile vide
	//partir d'une autre face non visisté
}*/


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
  	unsigned int nbVerts = 0;
	for (Vertex_iterator i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i) {
		++nbVerts;
	}
	std::cout << "Nombre de sommets: " << nbVerts << std::endl;
	
	unsigned int nbEdges = 0;
	for (Halfedge_iterator i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++i) {
		++nbEdges;
	}
	nbEdges /= 2;
	std::cout << "Nombre d'arêtes: " << nbEdges << std::endl;

	unsigned int nbFaces = 0;
	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i) {
		++nbFaces;
	}
	std::cout << "Nombre de faces: " << nbFaces << std::endl;
	unsigned int euler = nbVerts - nbEdges + nbFaces;
	unsigned int genus = (2 - euler) / 2;
	std::cout << "En supposant que le maillage contienne une unique surface sans bord, alors son genre est de " << genus << std::endl;
	
	float min, max, mean;
	std::map<Polyhedron::Facet_handle, double> perimetres = ComputePerimetre(mesh,min, max, mean);
	std::map<Polyhedron::Facet_handle, int> categories = ComputeTreshold(mesh, perimetres, mean);
	Draw(mesh, categories);
  
	return 0;
}
