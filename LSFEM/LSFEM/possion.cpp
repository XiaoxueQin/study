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

#include <AFEPack/DGFEMSpace.h>
#include <AFEPack/AMGSolver.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Operator.h>
#include <AFEPack/Functional.h>
#include <AFEPack/EasyMesh.h>
#include <AFEPack/HGeometry.h>

#include <lac/sparse_matrix.h>
#include <lac/sparsity_pattern.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>

#include <gsl/gsl_multifit.h>

#define PI (4.0*atan(1.0))
#define DIM 2
double u(const double *);
double f(const double *);

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

  

    /// 存放模版基函数 lambda, 第一维是单元在模版中编号, 第二维是模版编
    /// 号, 第三维是有限元基函数编号.
    std::vector<std::vector<gsl_vector* > > lambda(n_element);
    for(int i = 0; i < n_element; ++i)
	lambda[i].resize(ele_pattern[i].size());

    /// 准备模版基函数 lambda
    FEMSpace<double, DIM>::ElementIterator the_element = fem_space0.beginElement();
    FEMSpace<double, DIM>::ElementIterator end_element = fem_space0.endElement();


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
    }	

    // for ( int i = 0 ; i < n_element ; ++i)
    // 	for (int j = 0; j < ele_pattern[i].size(); ++j)
    // 	    for (int k = 0; k < 6; ++k)
		//	std::cout<<gsl_vector_get( lambda[i][j],k)<<std::endl;

   /// 先统计每一行非零元个数.
    std::vector<unsigned int> n_nz(n_element);
    // for (the_element = fem_space0.beginElement(); 
    // 	 the_element != end_element; 
    // 	 ++the_element) 
    // {
    // 	int ele_idx = the_element->index();
    // 	/// 每一单元的模版中任两个单元都会确定一个非零元.
    // 	int n_dof_pat = ele_pattern[ele_idx].size();
    // 	for (int i = 0; i < n_dof_pat; ++i)
    // 	{
    // 	    int row_id = ele_pattern[ele_idx][i]; 
    // 	    n_nz[row_id] += n_dof_pat;
    // 	    /// 防止小网格情形, 每行非零元个数超过行数.
    // 	    if (n_nz[row_id] > n_element)
    // 		n_nz[row_id] = n_element;
    // 	}
    // }    

 


// 统计非零元个数
    for (the_element = fem_space0.beginElement(); 
    	 the_element != end_element; 
    	 ++the_element) 
    {
    	int idx1=the_element->index();
    	GeometryBM &geo = the_element->geometry();
    	unsigned int n_boundary = geo.n_boundary();
	std::vector<int >idx2(n_boundary);
//	std::cout<<"OK:"<<n_boundary<<std::endl;
    	for( unsigned int i = 0; i < n_boundary; ++i )
    	{
    	    unsigned int dof_index = the_element->dof().at(0);
    	    unsigned int boundary_index = geo.boundary( i );
    	    DGElement<double, DIM> &dg_ele = fem_space0.dgElement( boundary_index );
    	    Element<double, DIM> *p_neigh = dg_ele.p_neighbourElement( 0 );
    	    if( p_neigh == &(*the_element) )
    		p_neigh = dg_ele.p_neighbourElement(1);
	    if( p_neigh == NULL )
		break;
	    else{
		idx1 = the_element->index();
		idx2[i] = p_neigh->index();
		int n_dof1 = ele_pattern[idx1].size();
		int n_dof2 = ele_pattern[idx2[i]].size();
		for (int j = 0; j < ele_pattern[idx1].size(); ++j)
		    n_nz[ele_pattern[idx1][j]] += n_dof2;
	    }
    	}
//	std::cout<<idx2[0]<<" , "<<idx2[1]<<" , "<<idx2[3]<<std::endl;
//		std::cout<<idx1<<"; "<<std::endl;	
    }

//非零元个数矫正
    for (int i = 0; i < n_element; ++i)
    {
	if (n_nz[i] > n_element)
	    n_nz[i]=n_element;
//	std::cout<<n_nz[i]<<std::endl;
    }

    SparsityPattern sp_mat(n_element, n_nz);

    // for (the_element = fem_space0.beginElement(); 
    // 	 the_element != end_element; 
    // 	 ++the_element) 
    // {
    // 	int ele_idx = the_element->index();
    // 	/// 每一单元的模版中任两个单元都会确定一个非零元.
    // 	int n_dof_pat = ele_pattern[ele_idx].size();
    // 	for (int i = 0; i < n_dof_pat; ++i)
    // 	{
    // 	    int row_id = ele_pattern[ele_idx][i];
    // 	    for (int j = 0; j < n_dof_pat; ++j)
    // 		sp_mat.add(row_id, ele_pattern[ele_idx][j]);
    // 	}
    // }

    //确定非零元在矩阵的位置，
    for (the_element = fem_space0.beginElement(); 
    	 the_element != end_element; 
    	 ++the_element) 
    {
    	int idx1=the_element->index();
    	GeometryBM &geo = the_element->geometry();
    	unsigned int n_boundary = geo.n_boundary();
	std::vector<int >idx2(n_boundary);
    	for( unsigned int i = 0; i < n_boundary; ++i )
    	{
    	    unsigned int dof_index = the_element->dof().at(0);
    	    unsigned int boundary_index = geo.boundary( i );
    	    DGElement<double, DIM> &dg_ele = fem_space0.dgElement( boundary_index );
    	    Element<double, DIM> *p_neigh = dg_ele.p_neighbourElement( 0 );
    	    if( p_neigh == &(*the_element) )
    		p_neigh = dg_ele.p_neighbourElement(1);
	    if( p_neigh == NULL )
		break;
	    else{
		idx1 = the_element->index();
		idx2[i] = p_neigh->index();
		for(int ii = 0; ii < ele_pattern[idx1].size(); ++ii)
		    for (int jj = 0; jj < ele_pattern[idx2[i]].size(); ++jj)
			sp_mat.add(ele_pattern[idx1][ii],ele_pattern[idx2[i]][jj]);
	    }
    	}
    }

    sp_mat.compress();



    SparseMatrix<double> stiff_matrix(sp_mat);
    //右端项
    Vector<double> right_hand_side(n_element);

    //矩阵第一部分（梯度内积）拼装/右端项拼装
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
			    cont = Jxw * innerProduct(basis_gradient[i][l], basis_gradient[j][l]) * gsl_vector_get(lambda[k][ii], i) * gsl_vector_get(lambda[k][jj],j);
			    row = ele_pattern[k][ii];
			    col = ele_pattern[k][jj];
			    stiff_matrix.add(row,col,cont);

			}
		    }
		    right_hand_side(row) += Jxw * basis_value[i][l] * f(q_point[l]) * gsl_vector_get(lambda[k][ii], i);  
		}
	    }
	}

//	std::cout<<stiff_matrix.diag_element(k)<<std::endl;
    }


    //矩阵第二部分（内部边界通量）拼装
    for (the_element = fem_space0.beginElement(); 
    	 the_element != end_element; 
    	 ++the_element) 
    {
    	int idx1=the_element->index();
    	GeometryBM &geo = the_element->geometry();
    	unsigned int n_boundary = geo.n_boundary();
	std::vector<int >idx2(n_boundary);
    	for( unsigned int i = 0; i < n_boundary; ++i )
	{
    	    unsigned int dof_index = the_element->dof().at(0);
    	    unsigned int boundary_index = geo.boundary( i );
    	    DGElement<double, DIM> &dg_ele = fem_space0.dgElement( boundary_index );
    	    Element<double, DIM> *p_neigh = dg_ele.p_neighbourElement( 0 );
    	    if( p_neigh == &(*the_element) )
    		p_neigh = dg_ele.p_neighbourElement(1);
	    if( p_neigh != NULL )
	    {
                 /// 计算边长
		double length = dg_ele.templateElement().volume();
		const QuadratureInfo<1>& dg_quad_info = dg_ele.findQuadratureInfo( 0 );
		std::vector<double> dg_jacobian = dg_ele.local_to_global_jacobian( dg_quad_info.quadraturePoint() );
		length *= dg_jacobian[ 0 ];
		/// 取外法向
		std::vector<Point<DIM> > dg_ele_ip = dg_ele.local_to_global( dg_quad_info.quadraturePoint() );
		std::vector<double> uno = unitOutNormal( dg_ele_ip[0], *the_element, dg_ele );
		std::vector<std::vector<double> > basis_value = the_element->basis_function_value(dg_ele_ip );

		//	std::cout<<basis_value[0][0]<<std::endl;
		std::cout<<gsl_vector_get( lambda[0][12],1)<<std::endl;
	    }
	
	}
    }
   
 

  
  

    //shifang hanshu kongjian
    for(int i = 0; i < n_element; ++i)
	for(int j = 0; j < ele_pattern[i].size(); ++j)
	    gsl_vector_free(lambda[i][j]);
 
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

