/**
 * @file   possion.cpp
 * @author Qin Xiaoxue 
 * @date   Thu Dec 17 13:44:21 2015
 * 
 * @brief  局部最小二乘求解微分方程模板.
 * 
 * 
 */

#include <iostream>
#include <fstream>
#include <list>

#include <AFEPack/AMGSolver.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Operator.h>
#include <AFEPack/Functional.h>
#include <AFEPack/EasyMesh.h>
#include <AFEPack/HGeometry.h>
#include <AFEPack/DGFEMSpace.h>

#include <lac/sparse_matrix.h>
#include <lac/sparsity_pattern.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>

#include <gsl/gsl_multifit.h>

#define PI (4.0*atan(1.0))
#define DIM 2

int main(int argc, char * argv[])
{
    /// 几何数.
    HGeometryTree<DIM> h_tree;
    h_tree.readEasyMesh(argv[1]);

    /// 两层网格链表, 0是数值解所在, 1是用于绘图.
    IrregularMesh<DIM> *irregular_mesh0;
    IrregularMesh<DIM> *irregular_mesh1;

    /// 生成数值解空间实际网格.
    irregular_mesh0 = new IrregularMesh<DIM>;
    irregular_mesh0->reinit(h_tree);
    irregular_mesh0->semiregularize();
    irregular_mesh0->regularize(false);
    RegularMesh<DIM> &mesh0 = irregular_mesh0->regularMesh();

    /// 生成绘图空间实际网格.
    irregular_mesh1 = new IrregularMesh<DIM>(*irregular_mesh0);
    irregular_mesh1->globalRefine(2);
    irregular_mesh1->semiregularize();
    irregular_mesh1->regularize(false);
    RegularMesh<DIM> &mesh1 = irregular_mesh1->regularMesh();

    /// 画出网格检测.
    mesh0.writeOpenDXData("D0.dx");
    mesh1.writeOpenDXData("D1.dx");

    /// 单元总数.
    int n_element = mesh0.n_geometry(DIM);
    /// 节点总数.
    int n_vertex = mesh0.n_geometry(0);

    /// 单元到顶点索引.
    std::vector<std::vector<int> > map_ele2vtx(n_element);

    /// 建立单元到顶点的索引.
    for (int i = 0; i < n_element; ++i)
    {
	GeometryBM &ele = mesh0.geometry(DIM, i); 
	int n_vtx = ele.n_vertex();
	map_ele2vtx[i].resize(n_vtx);
	for (int j = 0; j < n_vtx; ++j)
	    map_ele2vtx[i][j] = ele.vertex(j);
    }

    /// 顶点到单元的索引.
    std::vector<std::vector<int> > map_vtx2ele(n_vertex);
    /// 过渡变量, 统计每个顶点的相邻单元数.
    std::vector<int> n_vtx_per_ele(n_vertex, 0);
    for (int i = 0; i < n_element; ++i)
    {
	for (int j = 0; j < map_ele2vtx[i].size(); ++j)
	    n_vtx_per_ele[map_ele2vtx[i][j]]++;
    }
    /// 给每个顶点开辟正确的相邻单元总数空间.
    for (int i = 0; i < n_vertex; ++i)
    {
	map_vtx2ele[i].resize(n_vtx_per_ele[i]);
	/// 将其清零, 下面代码要利用此空间做统计.
	n_vtx_per_ele[i] = 0;
    }
    /// 然后赋入单元编号.
    for (int i = 0; i <  n_element; ++i)
    {
	for (int j = 0; j < map_ele2vtx[i].size(); ++j)
	{
	    /// 从单元读出顶点编号.
	    int idx_vtx = map_ele2vtx[i][j];
	    /// 将本单元编号放入对应顶点的单元索引中.
	    map_vtx2ele[idx_vtx][n_vtx_per_ele[idx_vtx]] = i;
	    /// 现在统计的是当前这个顶点下一次应该插入的单元索引的位置.
	    n_vtx_per_ele[idx_vtx]++;
	}
    }

    /// 单元模板.
    std::vector<std::vector<double> > ele_pattern(n_element);
    /// 构建单元模板.
    for (int i = 0; i < n_element; ++i)
    {

	int n_ele = 0;
	/// 通过一个双向链表, 动态以单元编号从小到大的次序插入.
	std::list<int> tmp_ele_pattern;
	int l = 0;
	for (int j = 0; j < map_ele2vtx[i].size(); ++j)
	{
	    for (int k = 0; k < map_vtx2ele[map_ele2vtx[i][j]].size(); ++k)
	    {
		int ele_idx2add = map_vtx2ele[map_ele2vtx[i][j]][k];
		std::list<int>::iterator pattern_iterator = tmp_ele_pattern.begin();
		std::list<int>::iterator pattern_end = tmp_ele_pattern.end();
		/// 判定是否有插入完成.
		bool isinserted = false;
		for (; pattern_iterator != pattern_end; ++pattern_iterator)
		    /// 若该单元已存在, 跳过.
		    if (*pattern_iterator == ele_idx2add)
		    {
			isinserted = true;
			break;
		    } /// 若当前单元编号大于等待插入单元编号, 则插入前一个位置.
		    else if (*pattern_iterator > ele_idx2add)
		    {
			isinserted = true;
			tmp_ele_pattern.insert(pattern_iterator, ele_idx2add);
			break;
		    }
		/// 之前未插入, 则插在最后.
		if (isinserted == false)
		{
		    tmp_ele_pattern.push_back(ele_idx2add);
		}
	    }
	    /// 将链表复制到数组索引.
	    ele_pattern[i].resize(tmp_ele_pattern.size());
	    std::list<int>::iterator pattern_iterator = tmp_ele_pattern.begin();
	    std::list<int>::iterator pattern_end = tmp_ele_pattern.end();
	    for (int j = 0; pattern_iterator != pattern_end; ++pattern_iterator, ++j)
		ele_pattern[i][j] = *pattern_iterator;
	}
	
    }
   
    ///建立含有基函数与高斯积分点的单元结构体
    TemplateGeometry<DIM> triangle_template_geometry;
    triangle_template_geometry.readData("triangle.tmp_geo");
    CoordTransform<DIM, DIM> triangle_coord_transform;
    triangle_coord_transform.readData("triangle.crd_trs");
    TemplateDOF<DIM> triangle_template_dof(triangle_template_geometry);
    triangle_template_dof.readData("triangle.2.tmp_dof");
    BasisFunctionAdmin<double, DIM, DIM> triangle_basis_function(triangle_template_dof);
    triangle_basis_function.readData("triangle.2.bas_fun");  

    UnitOutNormal<DIM> unit_out_normal;
    unit_out_normal.readData("triangle.out_nrm");

    std::vector<TemplateElement<double, DIM, DIM> > template_element(1);
    template_element[0].reinit(triangle_template_geometry,
			       triangle_template_dof,
			       triangle_coord_transform,
			       triangle_basis_function,
			       unit_out_normal);


    TemplateGeometry<DIM - 1> interval_template_geometry; /**< 区间单元（边）几何信息模板 */
    CoordTransform<DIM - 1, DIM> interval_to2d_coord_transform; /**< 1D 到 2D 的区间坐标变换模板 */
    std::vector<TemplateDGElement<DIM - 1, DIM> > dg_template_element; /**< 边界单元的 1D 参考单元模板 */

    DGFEMSpace<double, DIM> fem_space0;


    interval_template_geometry.readData("interval.tmp_geo");
    interval_to2d_coord_transform.readData("interval.to2d.crd_trs");
    dg_template_element.resize(1);
    dg_template_element[0].reinit(interval_template_geometry, interval_to2d_coord_transform);

    
    /// 有限元空间和网格的联系
    fem_space0.reinit(mesh0, template_element, dg_template_element);
    fem_space0.element().resize(n_element);
    for (int i = 0;i < n_element;i ++)
	fem_space0.element(i).reinit(fem_space0, i, 0);
    fem_space0.buildElement();
    fem_space0.buildDof();
    fem_space0.buildDofBoundaryMark();

    unsigned int n_side = mesh0.n_geometry(1);
    fem_space0.dgElement().resize(n_side);
    for( unsigned int i = 0; i < n_side; ++i )
         fem_space0.dgElement(i).reinit(fem_space0, i, 0);
    fem_space0.buildDGElement();

    int n_element1 = mesh1.n_geometry(2);
    
    /// 绘图有限元空间.
    FEMSpace<double, DIM> fem_space1(mesh1, template_element);	
    fem_space1.element().resize(n_element1);
    for (int i = 0;i < n_element1;i ++)
	fem_space1.element(i).reinit(fem_space1, i, 0);
    fem_space1.buildElement();
    fem_space1.buildDof();
    fem_space1.buildDofBoundaryMark();

    /// 暂时只有绘图数值解.
    FEMFunction<double, DIM> u_1(fem_space1);

    /// 存放模版基函数 lambda, 第一维是单元在模版中编号, 第二维是模版编
    /// 号, 第三维是有限元基函数编号.
    std::vector<std::vector<gsl_vector* > > lambda(n_element);
    for(int i = 0; i < n_element; ++i)
	lambda[i].resize(ele_pattern[i].size());

    /// 准备模版基函数 lambda
    FEMSpace<double, DIM>::ElementIterator the_element = fem_space0.beginElement();
    FEMSpace<double, DIM>::ElementIterator end_element = fem_space0.endElement();

    //画图函数
    std::vector<double > f_test(n_element);

    /// 构建在每一个单元上的所有基函数lambda
    for (;the_element != end_element; ++the_element) 
    {
        ///单元自由度
	int n_dof = the_element->n_dof();
        ///该单元的编号
	int ele_idx = the_element->index();
        ///该单元模板所含单元个数
	int n_ele_pat = ele_pattern[ele_idx].size();
	///模板的采样点
	std::vector<Point<DIM> > core_pnt(n_ele_pat);
	///利用0次精度的高斯积分点定义模板采样点
	for (int i = 0; i < n_ele_pat; ++i)
	{
	    Element<double, DIM> &pat_ele = fem_space0.element(ele_pattern[ele_idx][i]);
	    const QuadratureInfo<DIM>& quad_info = pat_ele.findQuadratureInfo(0);
	    core_pnt[i] = pat_ele.local_to_global(quad_info.quadraturePoint(0));
	}
	///计算出在此单元上的基函数   在该单元模板采样点的取值
	std::vector<std::vector<double> > basis_value = the_element->basis_function_value(core_pnt);
	///利用gsl 库中的最小二乘函数
	double xi, yi, ei, chisq;
	gsl_matrix *X, *cov;
	X = gsl_matrix_alloc (n_ele_pat, n_dof);
	cov = gsl_matrix_alloc (n_dof, n_dof);
	///超定方程组的定义
	for (int i = 0; i < n_ele_pat; ++i) 	    
	    for (int j = 0; j < n_dof; ++j)
		gsl_matrix_set (X, i, j, basis_value[j][i]);

	std::vector<gsl_vector* > yy(n_ele_pat);
	for(int i = 0 ; i < n_ele_pat; ++i)
	    yy[i]=gsl_vector_alloc(n_ele_pat);
  
	for(int j = 0 ; j < n_ele_pat; ++j)
	    lambda[ele_idx][j]=gsl_vector_alloc(n_dof);
       
	///定义不同函数在该单元下的采样值，并计算出系数
	for(int i = 0; i < n_ele_pat; ++i)
	{
	    for(int j = 0; j < n_ele_pat; ++j)
	    {
		if(j == i)
		    gsl_vector_set (yy[i], j, 1.0);
		else
		    gsl_vector_set (yy[i], j, 0.0);
	    }
	    gsl_multifit_linear_workspace *work 
		= gsl_multifit_linear_alloc (n_ele_pat, n_dof);
	    gsl_multifit_linear (X, yy[i], lambda[ele_idx][i], cov, &chisq, work);
	    gsl_multifit_linear_free (work);
	}
	gsl_matrix_free (X);
	for(int i = 0; i < n_ele_pat; ++i)
	    gsl_vector_free (yy[i]);
	gsl_matrix_free (cov);

        ///计算出所中心点的画图函数值
	std::vector<Point<DIM> > core_pnt1(n_element);
	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(0);
	core_pnt1[ele_idx] = the_element->local_to_global(quad_info.quadraturePoint(0));   
	f_test[ele_idx]=sin(PI*core_pnt1[ele_idx][0])*sin(PI*core_pnt1[ele_idx][1]);
    }	

    /// 先统计每一行非零元个数.
    std::vector<unsigned int> n_nz(n_element);
    for (the_element = fem_space0.beginElement(); 
	 the_element != end_element; 
	 ++the_element) 
    {
	int ele_idx = the_element->index();
	/// 每一单元的模版中任两个单元都会确定一个非零元.
	int n_dof_pat = ele_pattern[ele_idx].size();
	for (int i = 0; i < n_dof_pat; ++i)
	{
	    int row_id = ele_pattern[ele_idx][i]; 
	    n_nz[row_id] += n_dof_pat;
	    /// 防止小网格情形, 每行非零元个数超过行数.
	    if (n_nz[row_id] > n_element)
		n_nz[row_id] = n_element;
	}
    }

    SparsityPattern sp_mat(n_element, n_nz);

    for (the_element = fem_space0.beginElement(); 
	 the_element != end_element; 
	 ++the_element) 
    {
	int ele_idx = the_element->index();
	/// 每一单元的模版中任两个单元都会确定一个非零元.
	int n_dof_pat = ele_pattern[ele_idx].size();
	for (int i = 0; i < n_dof_pat; ++i)
	{
	    int row_id = ele_pattern[ele_idx][i];
	    for (int j = 0; j < n_dof_pat; ++j)
		sp_mat.add(row_id, ele_pattern[ele_idx][j]);
	}
    }
    
    sp_mat.compress();

    SparseMatrix<double> stiff_matrix(sp_mat);
    //右端项
    Vector<double> right_hand_side(n_element);


    for (the_element = fem_space0.beginElement(); 
	 the_element != end_element; ++the_element) 
    {
	int k = the_element->index();
	int n_ele_pat = ele_pattern[k].size();
	double volume = the_element->templateElement().volume();
	const QuadratureInfo<DIM>& quad_info = the_element->findQuadratureInfo(4);
	std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
	int n_quadrature_point = quad_info.n_quadraturePoint();
	std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
	std::vector<std::vector<std::vector<double> > > basis_gradient = the_element->basis_function_gradient(q_point);
	std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
	const std::vector<int>& element_dof = the_element->dof();
	int n_element_dof = the_element->n_dof();
	int row,col; 
	double cont;
	for (int l = 0; l < n_quadrature_point; ++l)
	{
	    double Jxw = quad_info.weight(l) * jacobian[l] * volume;
	  
	    for (int i = 0; i < n_element_dof; ++i)
	    {
		for (int ii = 0; ii < n_ele_pat; ++ii)
		{
		    for (int j = 0; j < n_element_dof; ++j)
		    {
			for (int jj = 0 ; jj < n_ele_pat; ++jj)
			{
			    cont = Jxw * innerProduct(basis_gradient[i][l], basis_gradient[j][l]) * gsl_vector_get(lambda[k][ii], l) * gsl_vector_get(lambda[k][jj],l);
			    //		    std::cout<<Jxw<<" , "<< innerProduct(basis_gradient[i][l], basis_gradient[j][l])<<" , "<<gsl_vector_get(lambda[k][ii], l)<<"  '"<< gsl_vector_get(lambda[k][jj],l)<<std::endl;
			    stiff_matrix.add(ele_pattern[k][ii],ele_pattern[k][jj],cont);

			}
		    }
		    //  right_hand_side(ele_pattern[k][ii]) += Jxw * basis_value[i][l] * f(q_point[l]) * gsl_vector_get(lambda[k][ii], l);  
		}
	    }
	}
	std::cout<<cont<<"   "<<std::endl;

	std::cout<<stiff_matrix.diag_element(k)<<std::endl;
    }

    

    std::ofstream output("solution.dat");

    output << "TITLE =\"LSFEM DATA\"" << std::endl;
    output << "VARIABLES = \"X\", \"Y\", \"u\"" << std::endl;
    output << "ZONE T=\"LSFEM\" NODES=" << n_element * 16 * 3 
	   << ", ELEMENTS=" << n_element * 16 
	   << ", DATAPACKING=POINT, ZONETYPE=FETRIANGLE" << std::endl;


    IrregularMesh<DIM, DIM>::RootIterator root_iterator = irregular_mesh1->beginRootElement();
    IrregularMesh<DIM, DIM>::RootIterator end_iterator = irregular_mesh1->endRootElement();

    //在细网格中画出函数值  
    for (; root_iterator != end_iterator; ++root_iterator)
    {
	int ele_idx = root_iterator->index;
	int n_chi = root_iterator->n_child;
	Element<double, DIM> the_element = fem_space0.element(ele_idx);
	int n_dof = the_element.n_dof();
	int n_ele_pat = ele_pattern[ele_idx].size();
        for (unsigned int i = 0; i < n_chi; ++i)
    	{ 
	    HElement<DIM, DIM> *chi = root_iterator->child[i]; 
	    int chi_ele_idx = chi->index;
	    int n_chi_chi = chi->n_child;
	    for (unsigned int j = 0; j < n_chi_chi; ++j)
	    {
		HElement<DIM, DIM> *chi_chi = chi->child[j];
		int chi_chi_ele_idx = chi_chi->index;
		Element<double, DIM> the_element_chi_chi = fem_space1.element(chi_chi_ele_idx);
		GeometryBM &geo = the_element_chi_chi.geometry();
		int n_vtx = geo.n_vertex();

		//粗网格的有限元基函数在细网格各个自由点的值	
		std::vector<std::vector<double> > value_test_pnt(n_vtx);
		for (int k = 0; k < n_vtx; ++k)
		    value_test_pnt[k] = the_element.basis_function_value(fem_space1.mesh().point(geo.vertex(k)));
		for (int k = 0; k < n_vtx; ++k)
		{
		    double z = 0;		    
                    
		    for(int ii = 0; ii < ele_pattern[ele_idx].size(); ++ii)
		     {
			 for (int l = 0; l < n_dof; ++l)
			     z += f_test[ele_pattern[ele_idx][ii]]* value_test_pnt[k][l]* gsl_vector_get(lambda[ele_idx][ii], l) ;
		     }
		    Point<DIM> &pnt = fem_space1.mesh().point(geo.vertex(k));
		    output << pnt[0] << "\t" << pnt[1] << "\t" << z << std::endl;
		}
	    }
	}
    }
    for (int i = 0; i < n_element * 16; ++i)
	output << i * 3 + 1 << " " << i * 3 + 2 << " " << i * 3 + 3 << std::endl;
    for(int i = 0; i < n_element; ++i)
	for(int j = 0; j < ele_pattern[i].size(); ++j)
	    gsl_vector_free(lambda[i][j]);

    output.close();
    return 0;
}

double u(const double * p)
{
    return sin(PI * p[0]) * sin(PI * p[1]);
};

double f(const double * p)
{
    return 2.0 * PI * PI * u(p);
};

