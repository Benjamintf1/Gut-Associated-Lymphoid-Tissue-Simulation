double generate_T(int x, int y, int , int height, double delta_space ){
	return 1000;
}

double generate_I(int x, int y, int width, int height, double delta_space ){

	return 0.0;
}


double generate_V(int x, int y, int width, int height, double delta_space ){

	if(x == width/2 && y == height/2){
		return 1500;
	}
	return 0.0;
}

double generate_birth(int x, int y, int width, int height, double delta_space ){

	return 0.00023;
}

int generate_tissues(int x, int y, int width, int height, double delta_space ){
	if(x < width/2) {
		return 0;
	} else {
		return 0;
	}
}


