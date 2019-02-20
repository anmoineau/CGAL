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
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;
typedef std::vector<double> Vect_double;
typedef std::vector<int> Vect_int;

double periMin = DBL_MAX;
double periMax = 0.0;
double periMean = 0.0;
double areaMin = DBL_MAX;
double areaMax = 0.0;

struct Color{
	float R;
	float G;
	float B;
	Color(float r, float g, float b): R(r), G(g), B(b){}
	Color(): R(0), G(0), B(0){}
};

void Draw(Polyhedron & P, std::map<Polyhedron::Facet_handle, int> categories){
	srand (time(NULL));
	std::map<int, Color> couleurs;
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
		Color c;
		if(couleurs.find(categories[i]) != couleurs.end()){
			c = couleurs[categories[i]];
		} else {
			couleurs[categories[i]] = Color(std::rand()/(float)RAND_MAX, std::rand()/(float)RAND_MAX, std::rand()/(float)RAND_MAX);
			c = couleurs[categories[i]];
		}
		file<<"  "<< c.R <<' '<<c.G<<' '<<c.B<<' '<<"0.75"<< std::endl;
        file << std::endl;
	}
}

std::map<Polyhedron::Facet_handle, double> ComputePerimetre(Polyhedron & P){
	
	double sum = 0;
	std::map<Polyhedron::Facet_handle, double> perimetres;
	for (Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		double perimetre = 0;
		Halfedge_facet_circulator j = i->facet_begin();
		CGAL_assertion( CGAL::circulator_size(j) >= 3);
		do {
            perimetre += sqrt(CGAL::squared_distance(j->vertex()->point(),j->opposite()->vertex()->point()));
        } while ( ++j != i->facet_begin());
		if(perimetre>periMax)
			periMax = perimetre;
		if(perimetre<periMin)
			periMin = perimetre;
		perimetres[i] = perimetre; 
		sum += perimetre;
	}
	periMean = sum/perimetres.size();
	std::cout << "Perimetre moyen: " << periMean << std::endl;
	return perimetres;
}

std::map<Polyhedron::Facet_handle, double> ComputeArea(Polyhedron &P) {
	std::map<Polyhedron::Facet_handle, double> aires;
	Vect_double dist;
	double p, aire;
	Facet_iterator f_it = P.facets_begin();
	while (f_it != P.facets_end()) {
		Halfedge_facet_circulator h_it = f_it->facet_begin();
		CGAL_assertion(CGAL::circulator_size(h_it) >= 3);
		do {
			dist.push_back(sqrt(CGAL::squared_distance(h_it->vertex()->point(), h_it->opposite()->vertex()->point())));
			++h_it;
		} while (h_it != f_it->facet_begin());
		p = (dist[0] + dist[1] + dist[2])/2;
		aire = sqrt(p * (p - dist[0]) * (p - dist[1]) * (p - dist[2]));
		aires[f_it] = aire;
		if (aire > areaMax) areaMax = aire;
		if (aire < areaMin) areaMin = aire;
		++f_it;
	}
	return aires;
}

std::map<Polyhedron::Facet_handle, int> ComputeTreshold(Polyhedron & P, std::map<Polyhedron::Facet_handle, double> & perimetres){
	std::map<Polyhedron::Facet_handle, int> categories;

    for (Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		categories[i] = (perimetres[i]> periMean)? 1 : 2;
	}
	return categories;
}

std::map<Polyhedron::Facet_handle, int> ComputeOtsu(Polyhedron &P, std::map<Polyhedron::Facet_handle, double> & perimetres ) {
	std::map<Polyhedron::Facet_handle, int> categories;
	Vect_int histo = Vect_int(64);
	Vect_int::iterator it = histo.begin();
	int section, threshold = 0;
	int nbFacets = perimetres.size();
	int sumA = 0, sumB = 0;
	int q1 = 0, q2;
	double u1, u2, current, max = 0.0;
	std::map<Polyhedron::Facet_handle, double> aire_map;
	aire_map = ComputeArea(P);
	while (it != histo.end()) {
		*it = 0;
		++it;
	}
	Facet_iterator f_it = P.facets_begin();
	while (f_it != P.facets_end()) {
		section = ((perimetres[f_it] - periMin) / (periMax- periMin)) * 64;
		histo[section] += ((aire_map[f_it] - areaMin) / (areaMax - areaMin)) * 100;
		++f_it;
	}
	for (int i = 0; i < 63; ++i) {
		sumA += i * histo[i];
	}
	for (int i = 0; i < 63; ++i) {
		q1 += histo[i];
		if (q1 == 0) continue;
		q2 = nbFacets - q1;
		sumB += i * histo[i];
		u1 = (double)sumB / (double)q1;
		u2 = (double)(sumA - sumB) / (double)q2;
		current = (double)q1 * (double)q2 * (u1 - u2) * (u1 - u2);
		if (current > max) {
			threshold = i;
			max = current;
		}
	}
	for (f_it = P.facets_begin() ; f_it != P.facets_end() ; ++f_it) {
		section = (perimetres[f_it] - periMin) / (periMax- periMin) * 64;
		categories[f_it] = (section > threshold)? 1 : 2;
	}
	return categories;
}

std::map<Polyhedron::Facet_handle, int> ModifyCategories(Polyhedron & P, std::map<Polyhedron::Facet_handle, int> & categories){
	int categorie;
	int subCategorie = 0;
	std::map<Polyhedron::Facet_handle, int> subCategories;
	std::map<Polyhedron::Facet_handle, bool> visite;
	std::queue<Polyhedron::Facet_handle> aVisiter;
	
	for (Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		visite[i] = false;
	}

	for (Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
		if(!visite[i]){
			subCategorie ++;
			aVisiter.push(i);
			visite[i]=true;
			categorie = categories[i];
			while(!aVisiter.empty()){
				Polyhedron::Facet_handle faceCour = aVisiter.front();
				aVisiter.pop();
				subCategories[faceCour]=subCategorie;
				Halfedge_facet_circulator j = faceCour->facet_begin();
				CGAL_assertion( CGAL::circulator_size(j) >= 3);
				do {
					if(!visite[j->opposite()->facet()] && categories[j->opposite()->facet()] == categorie){
						aVisiter.push(j->opposite()->facet());
						visite[j->opposite()->facet()] = true;
					}
				} while ( ++j != faceCour->facet_begin());
			}
		}	
	}
	return subCategories;
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
	
	std::map<Polyhedron::Facet_handle, double> perimetres = ComputePerimetre(mesh);
	std::map<Polyhedron::Facet_handle, int> categories = ComputeTreshold(mesh, perimetres);
	//std::map<Polyhedron::Facet_handle, int> categories = ComputeOtsu(mesh, perimetres);
	std::map<Polyhedron::Facet_handle, int> subCategories = ModifyCategories(mesh, categories);
	Draw(mesh, subCategories);
  
	return 0;
}
