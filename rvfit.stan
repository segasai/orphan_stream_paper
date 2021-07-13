// Author: Sergey Koposov ( Email : skoposov _AT_ ed.ac.uk)
// Changes
// July 2021 fix the of by one bug in the spline interpolation
// (does not affect the resutls of the paper)
functions{

	// get the vector of spacings between nodes	
	vector geths(int n_nodes, vector nodes)
	{
		int n = n_nodes -1;
		vector[n] hs;
		for (i in 1:n)
		{
			hs[i] = nodes[i+1] - nodes[i];
		}
		return hs;
	}

	// obtain the vector of spline coefficients given the location
	// of the nodes and values there 
	// We are using natural spline definition
	vector getcoeffs(int n_nodes, vector nodes, vector vals)
	{
		int n=n_nodes-1;
		vector[n] hi;
		vector[n] bi;
		vector[n-1] vi;
		vector[n-1] ui;
		vector[n_nodes] ret;
		vector[n-1] zs;
		matrix[n-1, n-1] M = rep_matrix(0, n-1, n-1);
		
		n = n_nodes-1;

		for (i in 1:n)
		{
			hi[i] = nodes[i+1]-nodes[i];
			bi[i] =  1/hi[i]*(vals[i+1]-vals[i]);
		}
		for (i in 2:n)
		{
			vi[i-1] = 2*(hi[i-1]+hi[i]);
			ui[i-1] = 6*(bi[i] - bi[i-1]);
		}
		for (i in 1:n-1)
		{
			M[i,i] = vi[i];
		}
		for (i in 1:n-2)
		{
			M[i+1,i] = hi[i+1];
			M[i,i+1] = hi[i+1];
		}
		//print (M)
		zs = M \ ui ; //mdivide_left_spd(M, ui);
		ret[1]=0;
		ret[n_nodes] =0;
		ret[2:n_nodes-1]=zs;
		
		return ret;
		
	}

	// Evaluate the spline, given nodes, values at the nodes
	// spline coefficients, locations of evaluation points 
	// and integer bin ids of each point	
	vector spline_eval(int n_nodes, vector nodes, 
	       vector vals, vector zs,
	       int n_dat, vector x, int[] i)
	{
		
		vector[n_nodes-1] h;
		vector[n_dat] ret;
		int i1[n_dat];
		for (ii in 1:n_dat)
		{
			i1[ii]=i[ii]+1;
		}
		h = geths(n_nodes, nodes);
		
		ret = ( 
		    zs[i1] ./ 6 ./ h[i] .* (x-nodes[i]).* (x-nodes[i]).*(x-nodes[i])+ 
		    zs[i]  ./ 6 ./ h[i] .* (nodes[i1]-x) .* (nodes[i1]-x) .* (nodes[i1]-x)+
		    (vals[i1] ./ h[i] - h[i] .* zs[i1] ./ 6).*(x-nodes[i])+
		    (vals[i] ./ h[i] - h[i] .* zs[i] ./ 6).*(nodes[i1]-x)
		    );
		return ret;
	}

	// find in which node interval we should place each point of the vector
	int[] findpos(int n_nodes, vector nodes, int n_dat, vector x)
	{
		int ret[n_dat];
		for (i in 1:n_dat)
		{
			for (j in 1:n_nodes-1)
			{
				if ((x[i]>=nodes[j]) && (x[i]<nodes[j+1]))
				{
					ret[i] = j;
				}
			}
		}
		return ret;
	}


}

data 
{
	// number of stars
	int n_stars;
	// number of spline nodes 
	int n_nodes;
	// number of spline nodes for the stream velocity dispersion
	int n_sig_nodes;
	// phi_1 of all the stream stars
	vector[n_stars] fi1;
	// radial velocities
	vector[n_stars] rv;
	// radial velocity errors
	vector[n_stars] erv;
	// nodes of the spline
	vector[n_nodes] fi1_nodes;
	// nodes of the spline for the velocity dispersion
	vector[n_sig_nodes] fi1_sig_nodes;
}
transformed data
{
	// the location of points in the spline nodes
	int node_ids[n_stars] = findpos(n_nodes, fi1_nodes, n_stars, fi1);
	int sig_node_ids[n_stars] = findpos(n_sig_nodes, fi1_sig_nodes, n_stars, fi1);
}

parameters
{
	// velocities at spline nodes
	vector[n_nodes] vels;
	// log of velocity dispersions 
	vector[n_sig_nodes] lsigs;

	// contamination fractions
	vector[n_nodes] lfracs; // these need to be transformed ex/(ex+1)
	// jacobian is also needed
	// background velocity mean
	real v0;
	// velocity dispersion of the contamination
	real<lower=0> sig0;
}

transformed parameters
{
}

model
{

	vector[n_stars] velmod;
	vector[n_stars] lsigmod;
	vector[n_stars] lfracmod;
	vector[n_stars] lfracmod1;
	vector[n_stars] lfracmod2;
	vector[n_stars] sigmod;
	vector[n_nodes] coeffs_vel= getcoeffs(n_nodes, fi1_nodes, vels);
	vector[n_nodes] coeffs_lfrac= getcoeffs(n_nodes, fi1_nodes, lfracs);
	//vector[n_nodes] coeffs_lsig = getcoeffs(n_nodes, fi1_nodes, lsigs);
	vector[n_sig_nodes] coeffs_lsig = getcoeffs(n_sig_nodes, fi1_sig_nodes, lsigs);
			
	// Priors
	sig0 ~ normal(110,20);
	v0 ~ normal (0,50);
	lsigs ~ normal(1.5,1.);	
	vels ~ normal(0,300);

	target += sum(lfracs-2*log1p_exp(lfracs)); // jacobian

	velmod  = spline_eval(n_nodes,fi1_nodes,
               vels, coeffs_vel,
	              n_stars, fi1, node_ids);
	lfracmod  = spline_eval(n_nodes, fi1_nodes,
               lfracs, coeffs_lfrac,
	              n_stars, fi1, node_ids);
	lsigmod  = spline_eval(n_sig_nodes, fi1_sig_nodes,
               lsigs, coeffs_lsig,
	              n_stars, fi1, sig_node_ids);
	
	sigmod = exp(lsigmod);
	lfracmod1 = lfracmod - log1p_exp(lfracmod);
	lfracmod2 = - log1p_exp(lfracmod);
	
	for (i in 1:n_stars)
	{
		target += log_sum_exp( lfracmod1[i] +
		       normal_lpdf(rv[i]|velmod[i], 
		       sqrt(sigmod[i] .*sigmod[i]+erv[i]*erv[i])),
	       lfracmod2[i] + 
	       	       normal_lpdf(rv[i]|v0,sqrt(sig0^2+erv[i]*erv[i])));
	}
}
