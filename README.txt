The following folder contain the main scripts, required functions, and a subset of the raw data used to model condensate-membrane interactions. All functions are self-contained in the one MATLAB script named "Deflection_Angle_Membrane_Interaction.m"

All code was written and tested on a MATLAB 2023/2024b Academic licence. Matlab signal analysis toolkit is required. When running the code - functions must be in the same directory as the main scripts. File directories to data will be system dependent and can be edited in the code - an example directory is provided.

Download time for required scripts, functions, sample data << 2 minutes from code repository. Scripts, functions, and sample data can be added to a MATLAB directory and run. 

Derivation of condensate membrane interaction physics and geometries are included in the supplemental modeling section.

The physical modeling of condensate-membrane interactions follows the general workflow described below:

	1. For a given condensate contact angle, solve for the above membrane and below membrane spherical cap 	radii 	assuming conservation of volume for different deflection angles: 

		i. doubleSphericalCapConstraint(delta_yDM_yDS)

		INPUTS:
		
		delta_yDM_yDS ranging between [-1 1] describes the intrinsic contact angle as a function of the adhesive work at the condensate-membrane interface using the Young's 			equation. Contact angles in radiant are obtained by taking arcos(delta_yDS_yDM).

		OUTPUTS:

		R_prime is an Nx1 vector containing the radius of the above membrane deflection radius for different above membrane deflection angles for a given intrinsic contact 			angle.

		phi is an Nx1 vector containing the different above membrane deflection angles. phi ranges between [0 intrinsic contact angle].

	2. For a given contact angle, calculate the free energy landscape and energy gradient for all possible geometries as determined by doubleSphericalCapConstraint.

		i. energyForGivenContactAngle(delta_yDM_yMS, k_bend, lambda)
		
		INPUTS: 

		delta_yDM_yDS ranging between [-1 1] describes the intrinsic contact angle as a function of the adhesive work at the condensate-membrane interface using the Young's 			equation. Contact angles in radiant are obtained by taking arcos(delta_yDS_yDM).

		k_bend is the non-dimensional membrane bending energy.

		lambda is the non-dimensional line tension.

		FUNCTION CALLS:

		doubleSphericalCapConstraint.

		surfaceEnergy(delta_yDM_yMS, yDS, lambda, theta, phi, R, k_bend. Calculates the free energy of condensate-membrane interactions provided a geometry and physical 			parameter. Refer to the supplemental modeling section for a description of all relevant physical terms. 

		OUTPUTS:
		
		E_surface is the free energy for a given contact angle for different above membrane deflection angles.
		
		theta_intrinsic, phi and R_prime are the geometric parameters as described above.

		ii. energyLandscape(theta_range, k_bend, lambda)
		
		populates a mesh-grid by iterating over energyForGivenContactAngle for different contact angles spaning [-pi pi].

		INPUTS:
		theta_range corresponds to a Nx1 vector of intrinsic contact angle values.
		k_bend is the non-dimensional membrane bending energy.
		lambda is the non-dimensional line tension.

		FUNCTION CALLS:
		energyForGivenContactAngle

		OUTPUTS:
		energyMatrix, phiMatrix, thetaMatrix correspond to NxM matrices containing normalize scalar free energies defined by their corresponding above membrane and intrinsic 			contact angles values. 

	3. Provided an energy landscape, we can numerically determine its gradient.

		i. gradientCountourPlotNoFitting(energyMatrix, phi_fraction, theta_normalized, delta_range, k_bend, lambda)

		INPUTS:
		energyMatrix is the scalar matrix containing the free energy values
		phi_fraction is the scalar matrix containing the above membrane deflection angles normalized by the intrinsic contact angle
		theta_normalized is the scalar matrix containing the intrinsic contact angle normalized by pi.
		delta_range is the vector contain all possible intrinsic contact angles.
		k_bend and lambda are physical parameters. 

		FUNCTION CALLS: 
		getCriticalPoints(delta_range, k_bend, lambda) determines the critical points for the free energy of a single contact angles as it deflects the membrane. 

		OUTPUTS:
		dE_Phi is the gradient of the free energy landscape with respect to phi (above membrane deflection angle) for a given intrinsic contact angle. We assume that condensate 		cannot transition to different contact angles.
	
	4. Provided experimental intrinsic and above membrane deflection angles, we can fit the gradient's stable solutions to determine estimates for k_bend and lambda

		i. experimentalSplineFit(expTheta, expPhi, delta_range, initial_guess)

		INPUTS:
		expTheta correspond to an Nx1 vector of measured intrinsic contact angles.
		expPhi correspond to an Nx1 vector of measured above membrane deflection angles.
		delta_range is the range of possible intrinsic contact angles.
		initial guess is a 2x1 vector containing the initial guess for k_bend and lambda.

		FUNCTION CALLS:
		splineResiduals(delta_range, expTheta, expPhi)
			fitSplineToCriticalPoints(delta_range, k_bend, lambda)
				getCriticalPoints(delta_range, k_bend, lambda)
		
		OUTPUTS:
		k_fit is the estimate for the non-dimensional membrane bending energy
		lambda_fit is the estimate for the non-dimensional lines tension

		ii. gradientCountourPlot(energyMatrix, phi_fraction, theta_normalized, delta_range, k_bend, lambda)

		functions the same as gradientContourPlotNoFitting, but also displays the stable solution corresponding to the k_fit and lambda_fit estimates determined in 				experimentalSplineFit.



				
	
	

	
