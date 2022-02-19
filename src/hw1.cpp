#include "hw1.h"
namespace algebra
{

    Matrix zeros(size_t n, size_t m)
    {	if (n==0||m==0)
		return Matrix{};
        Matrix a;
        {
            if (n < 0 || m < 0)
            {
                throw std::logic_error("logic error");
            }
            else
            {
                Matrix a(n, std::vector<double>(m, 0.0));
                return a;
            }
        }

        return a;
    }

    Matrix ones(size_t n, size_t m)
    {	if (n==0||m==0)
		return Matrix{};
        Matrix a;
        {
            if (n < 0 || m < 0)
            {
                throw std::logic_error("Eror");
            }
            else
            {
                Matrix a(n, std::vector<double>(m, 1.000));
                return a;
            }
        }

        return a;
    }

    Matrix random(size_t n, size_t m, double min, double max)
    {
        Matrix arr{zeros(n, m)};
        std::random_device rd;

        {
            if (max <= min)
            {
                throw std::logic_error("Eror");
            }
            else
            {
                for (size_t i{}; i < n; i++)
                {
                    for (size_t j{}; j < m; j++)
                    {
                        arr[i][j] = (static_cast<double>(rd())) / (rd.max() - rd.min()) * (max - min) + min;
                    }
                }
            }
        }

        return arr;
    }

    void show(const Matrix &matrix)
    {	if(matrix.empty())
		std::cout<<std::endl ;
	else{
        size_t n{(matrix).size()};
        size_t m{(matrix[0]).size()};
        for (size_t i{}; i < n; i++)
        {
            for (size_t j{}; j < m; j++)
                std::cout << std::setw(10) << std::setprecision(3) << matrix[i][j];
            std::cout << std::endl;
        }
	}
    }

    Matrix multiply(const Matrix &matrix, double c)
    {	
		if (matrix.empty())
            return Matrix{};
        size_t n{(matrix).size()};
        size_t m{(matrix[0]).size()};

        Matrix a{zeros(n, m)};
        for (size_t i{}; i < n; i++)
            for (size_t j{}; j < m; j++)
                a[i][j] = matrix[i][j] * c;

        return a;
    }

    Matrix multiply(const Matrix &matrix1, const Matrix &matrix2)
    {
        if (matrix1.empty() || matrix2.empty())
            return Matrix{};
        size_t n{(matrix1).size()};
        size_t m{(matrix2[0]).size()};
        size_t p{(matrix2).size()};
        double s{0.0};

        Matrix a{zeros(n, m)};
        {
            if ((p != matrix1[0].size()) || (matrix2.empty()) || (matrix1.empty()))
            {
                throw std::logic_error("Error");
            }
            else
            {
                for (size_t i{}; i < n; i++)
                    for (size_t j{}; j < m; j++)
                    {
                        for (size_t k{}; k < p; k++)
                        {
                            s += matrix1[i][k] * matrix2[k][j];
                        }
                        a[i][j] = s;
                        s = 0.0;
                    }
            }
        }

        return a;
    }
    Matrix sum(const Matrix &matrix, double c)
    {
        if (matrix.empty())
            return Matrix{};
        size_t n{(matrix).size()};
        size_t m{(matrix[0]).size()};
        Matrix a{zeros(n, m)};
        if (false)
        {
            throw std::logic_error("Error");
        }
        for (size_t i{}; i < n; i++)
            for (size_t j{}; j < m; j++)
                a[i][j] = matrix[i][j] + c;
        return a;
    }

    Matrix sum(const Matrix &matrix1, const Matrix &matrix2)
    {
        if (matrix1.empty() && matrix2.empty())
            return Matrix{};
        if (matrix1.empty() || matrix2.empty())
            throw std::logic_error("Error");
        size_t n{(matrix1).size()};
        size_t m{(matrix1[0]).size()};

        Matrix a{zeros(n, m)};
        {
            if ((matrix2[0].size() != m) || (matrix2.size() != n))
            {
                throw std::logic_error("Error");
            }
            else
            {
                for (size_t i{}; i < n; i++)
                    for (size_t j{}; j < m; j++)
                    {
                        a[i][j] = matrix1[i][j] + matrix2[i][j];
                    }
            }
        }

        return a;
    }

    Matrix transpose(const Matrix &matrix)
    {
        if (matrix.empty())
            return Matrix{};
        size_t n{(matrix).size()};
        size_t m{(matrix[0]).size()};
        Matrix a{zeros(m, n)};
        for (size_t i{}; i < n; i++)
            for (size_t j{}; j < m; j++)
                a[j][i] = matrix[i][j];
        return a;
    }

    Matrix minor(const Matrix &matrix, size_t n, size_t m)
    {	if (matrix.empty())
            return Matrix{};
        size_t row{(matrix).size()};
        size_t col{(matrix[0]).size()};
        Matrix a{zeros(row - 1, col - 1)};
        for (size_t i{}; i < row - 1; i++)
        {
            for (size_t j{}; j < col - 1; j++)
            {
                if (i < n && j < m)
                    a[i][j] = matrix[i][j];
                else if (i >= n && j < m)
                    a[i][j] = matrix[i + 1][j];
                else if (i < n && j >= m)
                    a[i][j] = matrix[i][j + 1];
                else
                    a[i][j] = matrix[i + 1][j + 1];
            }
        }

        return a;
    }

    double determinant(const Matrix &matrix)
    {
        if (matrix.empty())
            return 1;
        double s{0.000};
        int n{static_cast<int>((matrix).size())};
        if (n == 1)
            return matrix[0][0];
        {
            if (n != static_cast<int>(matrix[0].size()))
            {
                throw std::logic_error("Error");
            }
            else
            {
                if (n == 2)
                    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
                else
                {
                    for (int j{}; j < n; j++)
                    {
                        s = s + (static_cast<double>(-2 * (j % 2)) + 1) * matrix[0][j] * determinant(minor(matrix, 0, j));
                    }
                }
            }
        }

        return s;
    }

    Matrix inverse(const Matrix &matrix)
    {
        if (matrix.empty())
            return Matrix{};
        int n{static_cast<int>((matrix).size())};
        double det{determinant(matrix)};
        Matrix a{zeros(n, n)};
        {
            if (det == 0)
            {
                throw std::logic_error("Error");
            }
            else
            {
                for (int i{}; i < n; i++)
                {
                    for (int j{}; j < n; j++)
                        a[j][i] = (static_cast<double>(-2 * ((i + j) % 2) + 1)) * determinant(minor(matrix, i, j)) / det;
                }
            }
        }

        return a;
    }

    Matrix ero_swap(const Matrix &matrix, size_t r1, size_t r2)
    {	
		if (matrix.empty()){
			throw std::logic_error("Error");
							}
        size_t n{(matrix).size()};
        size_t m{(matrix[0]).size()};
        Matrix a{zeros(n, m)};

        {
            if (n - 1 < r1 || n - 1 < r2)
            {
                throw std::logic_error("Error");
            }
            else
            {
                for (size_t i{}; i < n; i++)
                    for (size_t j{}; j < m; j++)
                        if (i != r1 && i != r2)
                            a[i][j] = matrix[i][j];
                for (size_t j{}; j < m; j++)
                {
                    a[r1][j] = matrix[r2][j];
                    a[r2][j] = matrix[r1][j];
                }
            }
        }

        return a;
    }

    Matrix ero_multiply(const Matrix &matrix, size_t r, double c)
    {	
		if (matrix.empty()){
			throw std::logic_error("Error");
							}
        size_t n{(matrix).size()};
        size_t m{(matrix[0]).size()};
        Matrix a{zeros(n, m)};
        {
            if (n <= 0 || m <= 0)
            {
                throw std::logic_error("Error");
            }
            else
            {
                for (size_t i{}; i < n; i++)
                    for (size_t j{}; j < m; j++)
                        if (i != r)
                            a[i][j] = matrix[i][j];
                for (size_t j{}; j < m; j++)
                {
                    a[r][j] = c * matrix[r][j];
                }
            }
        }

        return a;
    }

    Matrix ero_sum(const Matrix &matrix, size_t r1, double c, size_t r2)
    {	
		if (matrix.empty()){
			throw std::logic_error("Error");
							}
        size_t n{(matrix).size()};
        size_t m{(matrix[0]).size()};
        Matrix a{zeros(n, m)};
        {
            if (n <= 0 || m <= 0)
            {
                throw std::logic_error("Error");
            }
            else
            {
                for (size_t i{}; i < n; i++)
                    for (size_t j{}; j < m; j++)
                        if (i != (r2))
                            a[i][j] = matrix[i][j];
                for (size_t j{}; j < m; j++)
                {
                    a[r2][j] = c * matrix[r1][j] + matrix[r2][j];
                }
            }
        }

        return a;
    }

    Matrix concatenate(const Matrix &matrix1, const Matrix &matrix2, int axis)
    {
        if (matrix1.empty() && matrix2.empty())
            return Matrix{};
        size_t n{(matrix1).size() + (axis ? 0 : (matrix2).size())};
        size_t m{(matrix1)[0].size() + ((!axis) ? 0 : (matrix2)[0].size())};
        Matrix a{zeros(n, m)};

        {
            if (((axis == 0) && ((matrix2)[0].size() != m)) || ((axis == 1) && ((matrix2).size() != n)))
            {
                throw std::logic_error("Error");
            }
            else
            {
                for (size_t i{}; i < (matrix1).size(); i++)
                    for (size_t j{}; j < (matrix1)[0].size(); j++)
                        a[i][j] = matrix1[i][j];
                if (axis == 0)
                {
                    for (size_t i{}; i < (matrix2).size(); i++)
                        for (size_t j{}; j < (matrix2)[0].size(); j++)
                            a[i + (matrix1).size()][j] = matrix2[i][j];
                }
                else
                {
                    for (size_t i{}; i < (matrix2).size(); i++)
                        for (size_t j{}; j < (matrix2)[0].size(); j++)
                            a[i][j + (matrix1[0]).size()] = matrix2[i][j];
                }
            }
        }

        return a;
    }

    Matrix upper_triangular(const Matrix &matrix)
    {
        if (matrix.empty())
            return Matrix{};
        size_t n{(matrix).size()};
        size_t m{(matrix[0]).size()};
        Matrix a{zeros(n, m)};
        if (n != m)
            throw std::logic_error("Error");
        for (size_t j{}; j < m; j++)
            for (size_t i{}; i < n; i++)
                a[i][j] = matrix[i][j];
        for (size_t j{}; j < m; j++)
        {
            //************************************************************ bounce question :
            if (a[j][j] == 0)
                for (size_t i{j + 1}; i < n; i++)
                    if (a[i][j] != 0)
                    {
                        a = ero_swap(a, i, j);
                        break;
                    }
					if (a[j][j] == 0) 
						throw std::logic_error("Eror");
            //************************************************************
            for (size_t i{j + 1}; i < n; i++)
            {
                a = ero_sum(a, j, -a[i][j] / a[j][j], i);
            }
        }
        return a;
    }
}