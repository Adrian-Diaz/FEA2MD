//Created by: Adrian Diaz
//Last edited by: Adrian Diaz
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#define PI 3.14159265
using namespace std;

int main()
{
	string line;
	int aiuc; //number of atoms in the unit cell
	int aiuc2;
	double cell[3][3]; //unit cell basis vectors and nanomesh basis vectors
	double cell2[3][3];
	double px, py, pz, pxn, pyn, pzn;
	double** atom_positions; //Will store the positions of the atoms in the unit cell
	int repeat[3]; //repeats of the unit cell in each dimension, and nanomesh unit cell scaling factor
	//compared to primitive unit cell, assumes orthogonal primitive unit cell
	int i, j, k, l, m, mm, mmm, atomn; //counter variables
	double unit_vector[2];
	int xflag,yflag,zflag;
	double box_delta=0;
  double box_size[3];
	double origin[3];
	double distance;
	double delr[3];
	double center_radius_sq;
	double wave_vector;
	double Amplitude = 0.01;
	double **nodes;
	int **element_sort, *element_type;
	double ***elements;
	int n_nodes, n_elements, dof_per_node, nodes_per_element,scale;

	/*************USER SETTINGS****************/
	dof_per_node = 3;
	nodes_per_element = 3;
	scale = 20;

	/***************END USER SETTINGS***********/
	ofstream myfile("Pt_cloak.txt"); //output filestream object for file output
	ofstream myfile2("atoms.xyz");
	ofstream myfile3("cellinfo.txt");
	ifstream infile("Basis-Pt.txt"); //input filestream object for file reading
	ifstream inmesh("Fullnew.msh"); //input filestream object for file reading
	//readin basis
	if (infile.is_open()) {
		infile >> line;
		cout << line << " ";
		infile >> line;
		cout << line << " ";
		infile >> line;
		cout << line << "\n";
		infile >> aiuc;
		aiuc2 = aiuc;
		cout << aiuc;
		for (i = 0; i<3; i++) {
			cout << "\n";
			for (j = 0; j<3; j++) {
				infile >> cell[i][j];
				cout << cell[i][j] << " ";
				cell2[i][j] = cell[i][j];
				//cell[i][j] = scale * cell[i][j];
			}
		}
		cout << "\n";

		infile >> repeat[0];

		infile >> repeat[1];

		infile >> repeat[2];


		cout << repeat[0] << " " << repeat[1] << " " << repeat[2] << "\n";
		atom_positions = new double*[aiuc2];
		for (i = 0; i<aiuc2; i++) {
			atom_positions[i] = new double[3];
		}
		for (i = 0; i<aiuc2; i++) {
			infile >> line;
			cout << line << " ";
			infile >> atom_positions[i][0];
			infile >> atom_positions[i][1];
			infile >> atom_positions[i][2];
			cout << atom_positions[i][0] << " " << atom_positions[i][1] << " " << atom_positions[i][2] << "\n";
		}

	}
	else {
		cout << "unable to open file";
	}
cout <<"reading FEA \n";
	//readin nodes
if (inmesh.is_open()) {
	//read in header
	  getline(inmesh,line);
		getline(inmesh,line);
		getline(inmesh,line);
		getline(inmesh,line);
		//read in number of nodes
		inmesh >> n_nodes;
		//allocate nodes array
		nodes = new double*[n_nodes];
		for (i = 0; i<n_nodes; i++) {
			nodes[i] = new double[dof_per_node];
		}
		cout << "\n";
		//read in nodal information
		for (i = 0; i<n_nodes; i++) {
			
			inmesh >> line;
			//cout << line << "\n";
			for (j = 0; j<dof_per_node; j++) {
				inmesh >> nodes[i][j];
				nodes[i][j]*=scale;
			}
			
		}
		//read in Endnodes and Elements strings
		getline(inmesh,line);
		//cout << line << "\n";
		getline(inmesh,line);
		//cout << line << "\n";
		getline(inmesh,line);
		cout << line << "\n";
		//read in number of elements
		inmesh >> n_elements;
		cout << n_elements << "\n";
		int index, dummy;
		//allocate element_sort array
		element_sort = new int*[n_elements];
		element_type = new int[n_elements];
		for (i = 0; i<n_elements; i++) {
			
			element_sort[i] = new int[nodes_per_element];
			//read element index
		  inmesh >>index;
			//cout << index << "\n";
			inmesh >>element_type[index-1];
			inmesh >>dummy;
			inmesh >>dummy;
			inmesh >>dummy;
			inmesh >> element_sort[index-1][0];
			inmesh >> element_sort[index-1][1];
			inmesh >> element_sort[index-1][2];
			
		}
			
	}
	else {
		cout << "unable to open file";
	}
 
 cout <<"bin setup \n";
 //fill per element data for convenience
 //allocate elements data
 elements = new double**[n_elements];
 //pointer assignment method
	/*
	for (i = 0; i<n_elements; i++)
		elements[i] = new double*[nodes_per_element];
  
  int index;
	
	for (i = 0; i<n_elements; i++) {
		for (j = 0; j>nodes_per_element; j++){
		index = element_sort[i][j]-1;
		elements[i][j] = nodes[index];
		}
	}
	*/

	//allocation method
	for (i = 0; i < n_elements; i++){
		elements[i] = new double*[nodes_per_element];
		for (j = 0; j < nodes_per_element; j++){
			elements[i][j] = new double[3];
		}
	}
  
  int index;
	for (i = 0; i<n_elements; i++) {
		for (j = 0; j < nodes_per_element; j++){
		index = element_sort[i][j]-1;
		elements[i][j][0] = nodes[index][0];
		elements[i][j][1] = nodes[index][1];
		elements[i][j][2] = nodes[index][2];
		}
	}

	//compute bounding box parameters
	double minx, miny, minz;
	double maxx, maxy, maxz;

	//initialize max and min values
	minx = maxx = nodes[0][0];
	miny = maxy = nodes[0][1];
	minz = maxz = nodes[0][2];
  
	//find min and max
	//loop over nodes and test coordinates for box box bounds
  for (i = 0; i<n_nodes; i++) {
	  if(nodes[i][0]>maxx) maxx = nodes[i][0];
		if(nodes[i][1]>maxy) maxy = nodes[i][1];
		if(nodes[i][2]>maxz) maxz = nodes[i][2];
		if(nodes[i][0]<minx) minx = nodes[i][0];
		if(nodes[i][1]<miny) miny = nodes[i][1];
		if(nodes[i][2]<minz) minz = nodes[i][2];
	}

  //set box size
	box_size[0] = maxx-minx;
	box_size[1] = maxy-miny;
	box_size[2] = maxz-minz;

  cout <<"bin allocate \n";
	//allocate bins
	int *bin_fill;
	int nbinx,nbiny,nbinz;
  nbinx = (int)(box_size[0]/cell[0][0])+1;
  nbiny = (int)(box_size[1]/cell[1][1])+1;
	nbinz = (int)(box_size[2]/cell[2][2])+1;
	bin_fill = new int[nbinx*nbiny*nbinz];
	cout <<"bin count: " << nbinx*nbiny*nbinz << "\n";

  //initialize bin_fill to zero
  cout <<"bin init \n";
	for(int zbin = 0; zbin < nbinx*nbiny*nbinz; zbin++){
		//cout << zbin << "\n";
		bin_fill[zbin] = 0;
	}
  
	int min_binx, min_biny, min_binz;
	int max_binx, max_biny, max_binz;
	int current_binx, current_biny, current_binz;
	//loop over elements and determine which bins they overlap
	cout <<"bin fill \n";
	for (i = 0; i<n_elements; i++) {
		//start by determining bounding box of bins this element may be concerned with
		//initialize bin limits
		cout << i << "\n";
		min_binx = max_binx = (int)((elements[i][0][0] - minx)/cell[0][0]);
		min_biny = max_biny = (int)((elements[i][0][1] - miny)/cell[1][1]);
		min_binz = max_binz = (int)((elements[i][0][2] - minz)/cell[2][2]);
   
		//test node for upper and lower bounds on bins to fill
		for (j = 0; j < nodes_per_element; j++) {
			current_binx = (int)((elements[i][j][0] - minx)/cell[0][0]);
			current_biny = (int)((elements[i][j][1] - miny)/cell[1][1]);
			current_binz = (int)((elements[i][j][2] - minz)/cell[2][2]);
			if(current_binx < min_binx) min_binx = current_binx;
			if(current_biny < min_biny) min_biny = current_biny;
			if(current_binz < min_binz) min_binz = current_binz;
			if(current_binx > max_binx) max_binx = current_binx;
			if(current_biny > max_biny) max_biny = current_biny;
			if(current_binz > max_binz) max_binz = current_binz;
		}
    int bin_index;
		/*
		if(min_binx!=0)
		min_binx-=1;
		if(min_biny!=0)
		min_biny-=1;
		if(min_binz!=0)
		min_binz-=1;
		if(max_binx!=nbinx-1)
		max_binx+=1;
		if(max_biny!=nbiny-1)
		max_biny+=1;
		if(max_binz!=nbinz-1)
		max_binz+=1;
    */

		if(element_type[i]==2){
		//define bounding surface (outward) normals for this element
    //define surface first then find normal
		double normals[3][2];
		double circuit_vectors[3][2];
		/*compute each normal using the standard cross product with the z unit vector as we traverse
		a counterclockwise circuit around the element*/
    int start_node, clockwise, node_select1, node_select2;
		double centroid[3], node2centroid[3][3], angles[3], norm[3];
		clockwise = 0;
    
		for(int idim=0; idim < 3; idim++){
		//compute element centroid
		centroid[idim] = (elements[i][0][idim] + elements[i][1][idim] + elements[i][2][idim])/3;
    
		//compute relative displacement of each node from the centroid
		node2centroid[0][idim] = elements[i][0][idim] - centroid[idim];
		node2centroid[1][idim] = elements[i][1][idim] - centroid[idim];
		node2centroid[2][idim] = elements[i][2][idim] - centroid[idim];
		}
    
		//angles subtended by nodes relative to centroid in current x-y coordinate system
		//return is from -pi/2 to pi/2
		norm[0] = sqrt(node2centroid[0][0]*node2centroid[0][0] + node2centroid[0][1]*node2centroid[0][1]);
		norm[1] = sqrt(node2centroid[1][0]*node2centroid[1][0] + node2centroid[1][1]*node2centroid[1][1]);
		norm[2] = sqrt(node2centroid[2][0]*node2centroid[2][0] + node2centroid[2][1]*node2centroid[2][1]);
		angles[0] = asin(node2centroid[0][1]/norm[0]);
		if(node2centroid[0][0]<0){
			if(angles[0]>0) angles[0] = PI - angles[0];
			else angles[0] = -PI - angles[0];
		}
		angles[1] = asin(node2centroid[1][1]/norm[1]);
		if(node2centroid[1][0]<0){
			if(angles[1]>0) angles[1] = PI - angles[1];
			else angles[1] = -PI - angles[1];
		}
		angles[2] = asin(node2centroid[2][1]/norm[2]);
		if(node2centroid[2][0]<0){
			if(angles[2]>0) angles[2] = PI - angles[2];
			else angles[2] = -PI - angles[2];
		}
		
		/*determine which node displacement traverses clockwise from node 1
		test for overlap of the periodic discontinuity in angular displacement by checking if 
    delta theta is greater than pi/2, the limit of angular displacement between two triangle vertices*/
		double theta12 = angles[1] - angles[0];
		double theta13 = angles[2] - angles[0];
		int discflag12, discflag13;
		discflag12 = discflag13 = 0;
		double pi = PI;
		if(theta12 >= pi || theta12 <= -pi) discflag12=1;
		if(theta13 >= pi || theta13 <= -pi) discflag13=1;

		if(theta12<0&&!discflag12) {
			node_select1 = 1;
			node_select2 = 2;
		}
		else if(theta13<0&&!discflag13) {
			node_select1 = 2;
			node_select2 = 1;
		}
		else if(discflag12) {
			node_select1 = 1;
			node_select2 = 2;
		}
    else {
			node_select1 = 2;
			node_select2 = 1;
		}

		//compute normal for the first part of the circuit using cross product
		for(int idim=0; idim < 2; idim++){
		circuit_vectors[0][idim] = elements[i][node_select1][idim] - elements[i][0][idim];
		circuit_vectors[1][idim] = elements[i][node_select2][idim] - elements[i][node_select1][idim];
		circuit_vectors[2][idim] = elements[i][0][idim] - elements[i][node_select2][idim];
		}
    normals[0][0] = -circuit_vectors[0][1];
		normals[0][1] = circuit_vectors[0][0];
		normals[1][0] = -circuit_vectors[1][1];
		normals[1][1] = circuit_vectors[1][0];
		normals[2][0] = -circuit_vectors[2][1];
		normals[2][1] = circuit_vectors[2][0];
		

		//loop through bounds on bins and test whether to fill more specifically for this element
		for(int ibin = min_binx; ibin <= max_binx; ibin++){
			for(int jbin = min_biny; jbin <= max_biny; jbin++){
				for(int kbin = min_binz; kbin <= max_binz; kbin++){
          bin_index = nbinx*nbiny*kbin + nbinx*jbin + ibin;

					//define triangular overlap calculation
					
					//cout << bin_index << "\n";
					  double bin_points[8][3];
						//define bin points w.r.t to origin in middle of bin
						
						//-x, -y, -z point
            bin_points[0][0] = minx + cell[0][0]*ibin;
						bin_points[0][1] = miny + cell[1][1]*jbin;
						bin_points[0][2] = minz + cell[2][2]*kbin;
            //x, -y, -z point
						bin_points[1][0] = minx + cell[0][0]*(ibin+1);
						bin_points[1][1] = miny + cell[1][1]*jbin;
						bin_points[1][2] = minz + cell[2][2]*kbin;
						//x, y, -z point
						bin_points[2][0] = minx + cell[0][0]*(ibin+1);
						bin_points[2][1] = miny + cell[1][1]*(jbin+1);
						bin_points[2][2] = minz + cell[2][2]*kbin;
						//-x, y, -z point
						bin_points[3][0] = minx + cell[0][0]*ibin;
						bin_points[3][1] = miny + cell[1][1]*(jbin+1);
						bin_points[3][2] = minz + cell[2][2]*kbin;
						//-x, -y, z point
            bin_points[4][0] = minx + cell[0][0]*ibin;
						bin_points[4][1] = miny + cell[1][1]*jbin;
						bin_points[4][2] = minz + cell[2][2]*(kbin+1);
            //x, -y, z point
						bin_points[5][0] = minx + cell[0][0]*(ibin+1);
						bin_points[5][1] = miny + cell[1][1]*jbin;
						bin_points[5][2] = minz + cell[2][2]*(kbin+1);
						//x, y, z point
						bin_points[6][0] = minx + cell[0][0]*(ibin+1);
						bin_points[6][1] = miny + cell[1][1]*(jbin+1);
						bin_points[6][2] = minz + cell[2][2]*(kbin+1);
						//-x, y, z point
						bin_points[7][0] = minx + cell[0][0]*ibin;
						bin_points[7][1] = miny + cell[1][1]*(jbin+1);
						bin_points[7][2] = minz + cell[2][2]*(kbin+1);

						/*projection of vector displacement, from bin point to a point on each surface, onto 
						the three normals must be negative for all three surfaces if the point lies inside the element. 
						Note that this may miss overlap of the element by a fraction of a unit cell; this is neglected.*/
						double displacement[3][3],inner_product;
						int sign_flag[3];
            for(int bintest=0; bintest < 8; bintest++){

              //displacement from point to the first element surface
						  for(int idim=0; idim < 3; idim++)
						  displacement[0][idim] = bin_points[bintest][idim]-elements[i][0][idim];
            
						  inner_product = normals[0][0]*displacement[0][0]+normals[0][1]*displacement[0][1];
						  if(inner_product < 0) sign_flag[0] = -1;
						  else sign_flag[0] = 1;

						  //displacement from point to the second element surface
						  for(int idim=0; idim < 3; idim++)
						  displacement[1][idim] = bin_points[bintest][idim]-elements[i][node_select1][idim];
            
						  inner_product = normals[1][0]*displacement[1][0]+normals[1][1]*displacement[1][1];
						  if(inner_product < 0) sign_flag[1] = -1;
						  else sign_flag[1] = 1;

						  //displacement from point to the third element surface
						  for(int idim=0; idim < 3; idim++)
						  displacement[2][idim] = bin_points[bintest][idim]-elements[i][0][idim];
            
						  inner_product = normals[2][0]*displacement[2][0]+normals[2][1]*displacement[2][1];
						  if(inner_product < 0) sign_flag[2] = -1;
						  else sign_flag[2] = 1;
            
						  if(sign_flag[0]==-1&&sign_flag[1]==-1&&sign_flag[2]==-1){
							  bin_fill[bin_index] = 1;
							  break;
						  }
						}
					}
				}
			}
		}
		else{
			for(int ibin = min_binx; ibin <= max_binx; ibin++){
			  for(int jbin = min_biny; jbin <= max_biny; jbin++){
				  for(int kbin = min_binz; kbin <= max_binz; kbin++){
						bin_index = nbinx*nbiny*kbin + nbinx*jbin + ibin;
						bin_fill[bin_index] = 1;
					}
				}
			}
		}
	}
  cout <<"end bin fill \n";
	//output MD model by testing bin positions
  //for(int zbin = 0; zbin < nbinx*nbiny*nbinz; zbin++){
		//cout << zbin << "\n";
		//if(bin_fill[zbin])
		//cout << "filled bin: " << zbin + 1 <<"\n";
	//}


	//define bound on number of necessary unit cell repeats to span the model
	repeat[0] = (int) box_size[0]/cell[0][0]+1;
	repeat[1] = (int) box_size[1]/cell[1][1]+1;
	repeat[2] = (int) box_size[2]/cell[2][2]+1;

	//output MD model
	cout <<"MD output \n";
	if (myfile.is_open()) {
		myfile << fixed << setprecision(8);
		myfile << "Platinum Cloak \n \n";
		myfile << repeat[0] * repeat[1] * repeat[2] << "        " << "atoms \n";
		myfile << "1 atom types \n \n";
		myfile << minx - box_delta<<"  " << maxx + box_delta<< "  xlo xhi \n";
		myfile << miny - box_delta<<"  " << maxy + box_delta<< "  ylo yhi \n";
		myfile << 0<<"  " << cell[2][2]<< "  zlo zhi \n";
		myfile << "Masses \n \n";
		myfile << "1  195.078 \n \n";
		myfile << "Atoms \n \n";
		origin[0] = (cell[0][0] * (repeat[0] - 0.5)-(cell[0][0] * -0.5))/2+ cell[0][0] * -0.5;
		origin[1] = (cell[1][1] * (repeat[1] - 0.5) - (cell[1][1] * -0.5)) / 2 + cell[1][1] * -0.5;
		origin[2] = (cell[2][2] * (repeat[2] - 0.5) - (cell[2][2] * -0.5)) / 2 + cell[2][2] * -0.5;
		cout <<"The origin " << origin[0] << " " << origin[1] << " " << origin[2] << "\n";

		atomn = 1;
	  int bin_index, binx, biny, binz;
		cout << repeat[0] << " " << repeat[1] << " " << repeat[2] << "\n";
					for (i = 0; i<repeat[0]; i++) {
						for (j = 0; j<repeat[1]; j++) {
							for (k = 0; k<repeat[2]; k++) {
                
								px = cell[0][0] * i + cell[1][0] * j + cell[2][0] * k + minx;
								py = cell[0][1] * i + cell[1][1] * j + cell[2][1] * k + miny;
								pz = cell[0][2] * i + cell[1][2] * j + cell[2][2] * k + minz;
                binx = (int)((px - minx)/cell[0][0]);
								biny = (int)((py - miny)/cell[1][1]);
								binz = (int)((pz - minz)/cell[2][2]);
				        bin_index = nbinx*nbiny*binz + nbinx*biny + binx;

								if (bin_fill[bin_index]) {
									//cout << "filled bin \n";
									for (l = 0; l < aiuc; l++) {
										px = cell[0][0] * i + cell[1][0] * j + cell[2][0] * k + atom_positions[l][0] + minx;
										py = cell[0][1] * i + cell[1][1] * j + cell[2][1] * k + atom_positions[l][1] + miny;
										pz = cell[0][2] * i + cell[1][2] * j + cell[2][2] * k + atom_positions[l][2] + minz;
										myfile <<std::left<<setw(10)<<atomn << setw(15)<< 1 <<setw(25)<<px<<setw(25)<<py<<setw(25)<<pz<<"\n";
										atomn = atomn + 1;
									}
								}
							}
						}
					}
				myfile.seekp(0, std::ios::beg);
				myfile << "Platinum Cloak \n \n";
		    myfile << atomn-1 << "        " << "atoms \n";
	}
	else {
		cout << "Unable to open file";
	}
	return 0;
}
