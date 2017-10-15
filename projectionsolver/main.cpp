#include <iostream>
#include <Eigen/Dense>

class NSsolver
{
public:
	NSsolver(int Nx, int Ny, double dx, double dy);
	enum BoundaryLocation { Up, Down, Left, Right };
	enum BoundaryType { Dirichlet, Neumann };
	void genMatU();
	void genLaplacianInnerMat(Eigen::MatrixXf & MatX);
	void genLaplacianBoundaryMat(Eigen::MatrixXf & MatX, enum BoundaryLocation bdloc, enum BoundaryType bdtp);
//	void genbU();

private:
	int Nx_;
	int Ny_;
	double dx_;
	double dy_;
	Eigen::VectorXd vecU_;
	Eigen::VectorXd vecV_;
	Eigen::VectorXd p_;
	Eigen::MatrixXf MatU_;
	Eigen::MatrixXf MatV_;
	Eigen::MatrixXf Matp_;
};

NSsolver::NSsolver(int Nx, int Ny, double dx, double dy)
	:Nx_(Nx),Ny_(Ny),dx_(dx),dy_(dy)
{
	vecU_.resize(Nx_*Ny_);
	vecV_.resize(Nx_*Ny_);
	p_.resize(Nx_*Ny_);
	MatU_.resize(Nx_*Ny_, Nx_*Ny_);
	MatV_.resize(Nx_*Ny_, Nx_*Ny_);
	Matp_.resize(Nx_*Ny_, Nx_*Ny_);
	for (int j = 0; j < Ny_; j++)
	{
		for (int i = 0; i < Nx_; i++)
		{
			int index = Nx_*j + i;
			vecU_(index) = 0;
			vecV_(index) = 0;
			p_(index) = 0;
		}
	}
	for (int j = 0; j < Nx_*Ny_; j++)
	{
		for (int i = 0; i < Nx_*Ny_; i++)
		{
			MatU_(i, j) = 0;
			MatV_(i, j) = 0;
			Matp_(i, j) = 0;
		}
	}
}

void NSsolver::genMatU()
{
	genLaplacianMat(MatU_);
	std::cout << MatU_;

}

void NSsolver::genLaplacianInnerMat(Eigen::MatrixXf & MatX)
{
	for (int j = 1; j < Ny_ - 1; j++)
	{
		for (int i = 1; i < Nx_ - 1; i++)
		{
			int index = Nx_*j + i;
			MatX(index, index) = MatX(index, index) - 2 / (dx_*dx_) - 2 / (dy_*dy_);
			MatX(index - 1, index) = MatX(index - 1, index) + 1 / (dx_*dx_);
			MatX(index + 1, index) = MatX(index + 1, index) + 1 / (dx_*dx_);
			MatX(index - Nx_, index) = MatX(index - Nx_, index) + 1 / (dy_*dy_);
			MatX(index + Nx_, index) = MatX(index + Nx_, index) + 1 / (dy_*dy_);
		}
	}
	for (int j = 1; j < Ny_ - 1; j++)
	{
		int i = 0;
		int index = Nx_*j + i;
		MatX(index, index) = MatX(index, index) - 2 / (dx_*dx_) - 2 / (dy_*dy_);
		MatX(index + 1, index) = MatX(index + 1, index) + 1 / (dx_*dx_);
		MatX(index - Nx_, index) = MatX(index - Nx_, index) + 1 / (dy_*dy_);
		MatX(index + Nx_, index) = MatX(index + Nx_, index) + 1 / (dy_*dy_);
	}
	for (int j = 1; j < Ny_ - 1; j++)
	{
		int i = Nx_ - 1;
		int index = Nx_*j + i;
		MatX(index, index) = MatX(index, index) - 2 / (dx_*dx_) - 2 / (dy_*dy_);
		MatX(index - 1, index) = MatX(index - 1, index) + 1 / (dx_*dx_);
		MatX(index - Nx_, index) = MatX(index - Nx_, index) + 1 / (dy_*dy_);
		MatX(index + Nx_, index) = MatX(index + Nx_, index) + 1 / (dy_*dy_);
	}
	for (int i = 1; i < Nx_ - 1; i++)
	{
		int j = 0;
		int index = Nx_*j + i;
		MatX(index, index) = MatX(index, index) - 2 / (dx_*dx_) - 2 / (dy_*dy_);
		MatX(index - 1, index) = MatX(index - 1, index) + 1 / (dx_*dx_);
		MatX(index + 1, index) = MatX(index + 1, index) + 1 / (dx_*dx_);
		MatX(index + Nx_, index) = MatX(index + Nx_, index) + 1 / (dy_*dy_);
	}
	for (int i = 1; i < Nx_ - 1; i++)
	{
		int j = Ny_ - 1;
		int index = Nx_*j + i;
		MatX(index, index) = MatX(index, index) - 2 / (dx_*dx_) - 2 / (dy_*dy_);
		MatX(index - 1, index) = MatX(index - 1, index) + 1 / (dx_*dx_);
		MatX(index + 1, index) = MatX(index + 1, index) + 1 / (dx_*dx_);
		MatX(index - Nx_, index) = MatX(index - Nx_, index) + 1 / (dy_*dy_);
	}
	{
		int i = 0;
		int j = 0;
		int index = Nx_*j + i;
		MatX(index, index) = MatX(index, index) - 2 / (dx_*dx_) - 2 / (dy_*dy_);
		MatX(index + 1, index) = MatX(index + 1, index) + 1 / (dx_*dx_);
		MatX(index + Nx_, index) = MatX(index + Nx_, index) + 1 / (dy_*dy_);
	}
	{
		int i = Nx_ - 1;
		int j = 0;
		int index = Nx_*j + i;
		MatX(index, index) = MatX(index, index) - 2 / (dx_*dx_) - 2 / (dy_*dy_);
		MatX(index - 1, index) = MatX(index - 1, index) + 1 / (dx_*dx_);
		MatX(index + Nx_, index) = MatX(index + Nx_, index) + 1 / (dy_*dy_);
	}
	{
		int i = 0;
		int j = Ny_ - 1;
		int index = Nx_*j + i;
		MatX(index, index) = MatX(index, index) - 2 / (dx_*dx_) - 2 / (dy_*dy_);
		MatX(index + 1, index) = MatX(index + 1, index) + 1 / (dx_*dx_);
		MatX(index - Nx_, index) = MatX(index - Nx_, index) + 1 / (dy_*dy_);
	}
	{
		int i = Nx_ - 1;
		int j = Ny_ - 1;
		int index = Nx_*j + i;
		MatX(index, index) = MatX(index, index) - 2 / (dx_*dx_) - 2 / (dy_*dy_);
		MatX(index - 1, index) = MatX(index - 1, index) + 1 / (dx_*dx_);
		MatX(index - Nx_, index) = MatX(index - Nx_, index) + 1 / (dy_*dy_);
	}
}

void NSsolver::genLaplacianBoundaryMat(Eigen::MatrixXf & MatX, enum BoundaryLocation bdloc, enum BoundaryType bdtp)
{
	if (bdloc == Up)
	{

	}
}




int main()
{
	const int M = 3;
	const int dx = 1;
	const int N = 3;
	const int dy = 1;
	double U[M][N];
	double V[M][N];
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < M; i++)
		{
		}
	}
	NSsolver nstest(M, N, dx, dy);
	nstest.genMatU();
	std::cin.get();
}