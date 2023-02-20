#include <TinyAD/ScalarFunction.hh>
#include <iostream>
#include <igl/matlab_format.h>
int main()
{
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> V(3,2);
  V<<0,1,2,3,4,5;
  auto func = TinyAD::scalar_function<2>(TinyAD::range(V.rows()));
  func.add_elements<1>(TinyAD::range(V.rows()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element)
  {
    // Evaluate element using either double or TinyAD::Double
    using T = TINYAD_SCALAR_TYPE(element);
    // Get variable 2D vertex positions
    Eigen::Vector2<T> v = element.variables(element.handle);
    return 0.5*v(0)*v(0) + 10*0.5*v(1)*v(1);
  });
  Eigen::VectorXd x = func.x_from_data([&] (int i) { return V.row(i); });
  auto [f, g, H] = func.eval_with_derivatives(x);
  std::cout<<igl::matlab_format(V,"V")<<std::endl;
  std::cout<<igl::matlab_format(x,"x")<<std::endl;
  std::cout<<igl::matlab_format(g,"g")<<std::endl;
  std::cout<<igl::matlab_format(Eigen::MatrixXd(H),"H")<<std::endl;
}
