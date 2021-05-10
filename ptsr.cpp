#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>

#include <CGAL/Vector_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/centroid.h>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <algorithm>
#include <queue>
#include <iostream>
#include <sys/time.h>
#include <boost/bind.hpp>
#include <CGAL/Object.h>
#include "nanoflann.hpp"
//#include "avltree.hpp"



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Simple_cartesian<CGAL::Quotient<CGAL::MP_Float> >  EK; // used for kernal conversion to be used for CGAL::intersection

typedef CGAL::Triangulation_cell_base_with_info_3<int, K> Cb;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb,Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<K>   Delaunay;
typedef CGAL::Triangle_3<K>                 Triangle;
typedef CGAL::Point_3<K>                    Point_3;
//typedef CGAL::Tetrahedron_3<K>              Tetrahedron;
typedef CGAL::Segment_3<K>                  Segment;
typedef CGAL::Line_2<K>                     Line2;
typedef CGAL::Segment_2<K>                  Segment2;
typedef CGAL::Triangle_2<K>                 Triangle2;
typedef CGAL::Point_2<K>                    Point2;
typedef Delaunay::Point                     Point;
typedef Delaunay::Edge                      Edge;
typedef Delaunay::Facet                     Facet;
typedef K::Plane_3                          Plane;
typedef CGAL::Cartesian_converter<K,EK>                         K_to_EK;
typedef CGAL::Cartesian_converter<EK,K>                         EK_to_K;

typedef unsigned long long timestamp_t;
using namespace nanoflann;

Delaunay T;


//The struct PointCloud is used to create the kd-tree for the point cloud.

struct PointCloud
{
    struct Point
    {
        float  x, y, z;
    };

    std::vector<Point>  pts;

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const { return pts.size(); }

    inline float kdtree_distance(const float *p1, const size_t idx_p2, size_t /*size*/) const
    {
        const float d0 = p1[0] - pts[idx_p2].x;
        const float d1 = p1[1] - pts[idx_p2].y;
        const float d2 = p1[2] - pts[idx_p2].z;
        return d0*d0 + d1*d1 + d2*d2;
    }

    inline float kdtree_get_pt(const size_t idx, int dim) const
    {
        if (dim == 0) return pts[idx].x;
        else if (dim == 1) return pts[idx].y;
              else return pts[idx].z;
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }

  };

typedef KDTreeSingleIndexAdaptor< L2_Simple_Adaptor< float, PointCloud >, PointCloud, 3 >  KDTree;

std::vector<Point> vertice; // to hold the points
std::vector<Point> vert_for_plane_fit;
 std::vector<Triangle> original_faces;
 typedef std::vector<Triangle> FaceList;

//const size_t num_nbr;// = 24;
size_t num_nbr;
const size_t vert_indx = 1;
 class Nbr_Tri_of_Point
 {
public:
   int done;
   FaceList facelist;
   FaceList black_list;
 };

struct Del_Faces_Detail
{
    Triangle tri;
    double circum_rad;
};

 Del_Faces_Detail del_face_detail[5000];

std::vector<Triangle> mesh;

 int check_overlap_trn(Triangle t1,Triangle t2, Plane plane)
 {

	    K_to_EK to_exact;
	    EK_to_K to_inexact;

	    EK::Point_2 pp;
 	    EK::Segment_2 seg;
 	    EK::Triangle_2 tr;

      std::vector<EK::Point_2> poly;

	    Point2 p_ref = plane.to_2d(t1.vertex(0));
     	Point2 p1_1 = plane.to_2d(t1.vertex(1));
      Point2 p1_2 = plane.to_2d(t1.vertex(2));

      Point2 p2_1 = plane.to_2d(t2.vertex(1));
      Point2 p2_2 = plane.to_2d(t2.vertex(2));
      Triangle2 tr1 = Triangle2(p_ref,p1_1,p1_2);
      Triangle2 tr2 = Triangle2(p_ref,p2_1,p2_2);
      EK::Triangle_2 tr11=to_exact(tr1);
      EK::Triangle_2 tr22=to_exact(tr2);


      CGAL::cpp11::result_of<EK::Intersect_2(EK::Triangle_2, EK::Triangle_2)>::type intr = CGAL::intersection(tr11,tr22);

      if(CGAL::assign(pp,intr) || CGAL::assign(seg,intr))
	     {
   	        return 0;
	     }
	    else if(CGAL::assign(tr,intr) || CGAL::assign(poly,intr))
	         {
		             return 1;  //triangle or polygon
	         }
      // return 0;
    }


 int test(Point ref,std::vector<Triangle> face_list,Point p1,Point p2)
 {

      if(face_list.size() == 0)
          return 0;
      std::vector<Point> plane_fit_ver;
      Point p_proj1,p_proj2,p_ref,temp_p_proj1,temp_p_proj2,temp_p_ref_1;
      int intrsect;
      Plane plane;
      for(int l=0;l<face_list.size();l++)
      {
          plane_fit_ver.push_back(face_list[l].vertex(0));
          plane_fit_ver.push_back(face_list[l].vertex(1));
          plane_fit_ver.push_back(face_list[l].vertex(2));
      }
      plane_fit_ver.push_back(p1);
      plane_fit_ver.push_back(p2);
      //plane_fit_ver.push_back(tr.vertex(0));

      linear_least_squares_fitting_3(plane_fit_ver.begin(),plane_fit_ver.end(),plane,CGAL::Dimension_tag<0>());

      p_proj1= plane.projection(p1);
      p_proj2 = plane.projection(p2);
      p_ref = plane.projection(ref);


      for(int j=0;j<face_list.size();j++)
      {
           if(ref == face_list[j].vertex(0))
           {
             temp_p_proj1= plane.projection(face_list[j].vertex(1));
             temp_p_proj2 = plane.projection(face_list[j].vertex(2));
           }
           else if(ref == face_list[j].vertex(1))
           {
             temp_p_proj1= plane.projection(face_list[j].vertex(0));
             temp_p_proj2 = plane.projection(face_list[j].vertex(2));
           }
           else
           {
             temp_p_proj1= plane.projection(face_list[j].vertex(0));
             temp_p_proj2 = plane.projection(face_list[j].vertex(1));
           }


           intrsect = check_overlap_trn(Triangle(p_ref,temp_p_proj1,temp_p_proj2),Triangle(p_ref,p_proj1,p_proj2),plane);

           if(intrsect)
           {
               return 1;
           }

      }

      return 0;

 }

 int check_manifold(Point ref, std::vector<Triangle> original_faces)
 {
       if(original_faces.size()==0)
         return 0;
        int charge = 0,face_size,mfld =0;
       Point p_ref,p_proj1,p_proj2;
       double x,y,z;
       face_size = original_faces.size();
       std::vector<Point> mid_edge;

       for(int j=0;j<face_size;j++)
       {
           mid_edge.push_back(original_faces[j].vertex(1));
           mid_edge.push_back(original_faces[j].vertex(2));
       }

       for(int j=0;j<mid_edge.size();j++)
     {     mfld = 0; //
        charge = charge + 1;
        for(int i=0;i<mid_edge.size();i++)
            if(mid_edge[j] == mid_edge[i] && i!=j)
            {
                charge = charge - 1;
                mfld++;
             }
        if(mfld >= 2)
           break;
      }

     if(mfld<2)
         return charge;
     else
          return -13;
 }


 bool triEqual(Triangle t1, Triangle t2)
 {
 	Point_3 c10 = t1.vertex(0);
 	Point_3 c11 = t1.vertex(1);
 	Point_3 c12 = t1.vertex(2);

 	Point_3 c20 = t2.vertex(0);
 	Point_3 c21 = t2.vertex(1);
 	Point_3 c22 = t2.vertex(2);

 	int number = 0;
 	if (c10 == c20 || c10 == c21 || c10 == c22)
 		number++;
 	if (c11 == c20 || c11 == c21 || c11 == c22)
 		number++;
 	if (c12 == c20 || c12 == c21 || c12 == c22)
 		number++;

 	if (number == 3)
 		return true;
 	return false;
 }

 std::vector<Triangle> recon(Point ref,Del_Faces_Detail *del_face_detail,int counter,std::vector<Triangle> fixed_faces, KDTree* pcdIndex,Nbr_Tri_of_Point* nbr_tri_of_point,int sgn,Nbr_Tri_of_Point* nbr_tri_buffer)
 {

        Plane plane;
        int flg = 1,test_1 = -1,test_2 =-1,test_3 =-1, test_4 = -1;
        int intrsect = 0;
        Point p_ref,p_proj1,p_proj2,p1,p2,temp_p_proj1,temp_p_proj2,t1,t2;

        std::vector<size_t> ret_index(num_nbr);
        std::vector<size_t> ret_index_1(num_nbr);
        std::vector<size_t> ret_index_2(num_nbr);
        std::vector<float> out_dist_sqr(num_nbr);

        std::vector<Triangle> buffer;
        std::vector<Triangle> additional_faces;
        vert_for_plane_fit.clear();
        original_faces.clear();



        vert_for_plane_fit.push_back(ref);
        int cros1 =0,cros2 =0,cros3 =0;
        int t11,t22,t33;
        std::vector<Triangle> tri_align_vert;
        std::vector<Triangle> temp_tri;
        if(fixed_faces.size() == 0)
        {
              int tr = 1,id=0;
              while(tr && id<counter)
              {
                  float queryPoint[3] = {del_face_detail[id].tri.vertex(0).x(), del_face_detail[id].tri.vertex(0).y(),del_face_detail[id].tri.vertex(0).z()};
                  pcdIndex->knnSearch(&queryPoint[0], vert_indx, &ret_index[0], &out_dist_sqr[0]);
                  float queryPoint1[3] = {del_face_detail[id].tri.vertex(1).x(), del_face_detail[id].tri.vertex(1).y(),del_face_detail[id].tri.vertex(1).z()};
                  pcdIndex->knnSearch(&queryPoint1[0], vert_indx, &ret_index_1[0], &out_dist_sqr[0]);
                  float queryPoint2[3] = {del_face_detail[id].tri.vertex(2).x(), del_face_detail[id].tri.vertex(2).y(),del_face_detail[id].tri.vertex(2).z() };
                  pcdIndex->knnSearch(&queryPoint2[0], vert_indx, &ret_index_2[0], &out_dist_sqr[0]);
                  cros1 =0;cros2 =0;cros3 =0;
                  cros1 = test(del_face_detail[id].tri.vertex(0),nbr_tri_of_point[ret_index[0]].facelist,del_face_detail[id].tri.vertex(1),del_face_detail[id].tri.vertex(2));
                  tri_align_vert.clear();
                  tri_align_vert = nbr_tri_of_point[ret_index[0]].facelist;
                  tri_align_vert.push_back(Triangle(del_face_detail[id].tri.vertex(0),del_face_detail[id].tri.vertex(1),del_face_detail[id].tri.vertex(2)));
                  t11 = check_manifold(del_face_detail[id].tri.vertex(0),tri_align_vert);

                  cros2 = test(del_face_detail[id].tri.vertex(1),nbr_tri_of_point[ret_index_1[0]].facelist,del_face_detail[id].tri.vertex(2),del_face_detail[id].tri.vertex(0));
                  tri_align_vert.clear();
                  tri_align_vert = nbr_tri_of_point[ret_index_1[0]].facelist;
                  tri_align_vert.push_back(Triangle(del_face_detail[id].tri.vertex(1),del_face_detail[id].tri.vertex(2),del_face_detail[id].tri.vertex(0)));
                  t22 = check_manifold(del_face_detail[id].tri.vertex(1),tri_align_vert);

                  cros3 = test(del_face_detail[id].tri.vertex(2),nbr_tri_of_point[ret_index_2[0]].facelist,del_face_detail[id].tri.vertex(0),del_face_detail[id].tri.vertex(1));
                  tri_align_vert.clear();
                  tri_align_vert = nbr_tri_of_point[ret_index_2[0]].facelist;
                  tri_align_vert.push_back(Triangle(del_face_detail[id].tri.vertex(2),del_face_detail[id].tri.vertex(0),del_face_detail[id].tri.vertex(1)));
                  t33 = check_manifold(del_face_detail[id].tri.vertex(2),tri_align_vert);


                 if(cros1 || cros2 || cros3 || nbr_tri_of_point[ret_index[0]].done ==13 || nbr_tri_of_point[ret_index_1[0]].done ==13 || nbr_tri_of_point[ret_index_2[0]].done ==13 )
                      id++;
                 else if(t11 < 0 || t22 < 0 || t33 < 0)
                      id++;
                  else
                      tr=0;
                }

              if(id<counter)
              {
                  if(ref == del_face_detail[id].tri.vertex(0))
                  {
                      vert_for_plane_fit.push_back(Point(del_face_detail[id].tri.vertex(1)));
                      vert_for_plane_fit.push_back(Point(del_face_detail[id].tri.vertex(2)));
                      nbr_tri_buffer[ret_index_1[0]].facelist.push_back(Triangle(del_face_detail[id].tri.vertex(1),del_face_detail[id].tri.vertex(2),del_face_detail[id].tri.vertex(0)));
                      nbr_tri_buffer[ret_index_2[0]].facelist.push_back(Triangle(del_face_detail[id].tri.vertex(2),del_face_detail[id].tri.vertex(0),del_face_detail[id].tri.vertex(1)));
                  }
                  else if(ref == del_face_detail[id].tri.vertex(1))
                  {
                      vert_for_plane_fit.push_back(Point(del_face_detail[id].tri.vertex(0)));
                      vert_for_plane_fit.push_back(Point(del_face_detail[id].tri.vertex(2)));
                      nbr_tri_buffer[ret_index[0]].facelist.push_back(Triangle(del_face_detail[id].tri.vertex(0),del_face_detail[id].tri.vertex(1),del_face_detail[id].tri.vertex(2)));
                      nbr_tri_buffer[ret_index_2[0]].facelist.push_back(Triangle(del_face_detail[id].tri.vertex(2),del_face_detail[id].tri.vertex(0),del_face_detail[id].tri.vertex(1)));
                  }
                  else
                  {
                      vert_for_plane_fit.push_back(Point(del_face_detail[id].tri.vertex(1)));
                      vert_for_plane_fit.push_back(Point(del_face_detail[id].tri.vertex(0)));
                      nbr_tri_buffer[ret_index_1[0]].facelist.push_back(Triangle(del_face_detail[id].tri.vertex(1),del_face_detail[id].tri.vertex(2),del_face_detail[id].tri.vertex(0)));
                      nbr_tri_buffer[ret_index[0]].facelist.push_back(Triangle(del_face_detail[id].tri.vertex(0),del_face_detail[id].tri.vertex(2),del_face_detail[id].tri.vertex(1)));
                  }

                  original_faces.push_back(Triangle(ref,vert_for_plane_fit[1],vert_for_plane_fit[2]));
                  additional_faces.push_back(Triangle(ref,vert_for_plane_fit[1],vert_for_plane_fit[2]));
             }
             else
                return additional_faces;
          }
          else if(fixed_faces.size() > 0)
          {
              for(int i = 0; i<fixed_faces.size();i++)
              {
                  if(ref == fixed_faces[i].vertex(0))
                  {
                      vert_for_plane_fit.push_back(Point(fixed_faces[i].vertex(1)));
                      vert_for_plane_fit.push_back(Point(fixed_faces[i].vertex(2)));
                      original_faces.push_back(Triangle(ref,Point(fixed_faces[i].vertex(1)),Point(fixed_faces[i].vertex(2))));

                   }
                   else if(ref == fixed_faces[i].vertex(1))
                   {
                          vert_for_plane_fit.push_back(Point(fixed_faces[i].vertex(0)));
                          vert_for_plane_fit.push_back(Point(fixed_faces[i].vertex(2)));
                          original_faces.push_back(Triangle(ref,Point(fixed_faces[i].vertex(0)),Point(fixed_faces[i].vertex(2))));
                   }
                  else
                   {
                          vert_for_plane_fit.push_back(Point(fixed_faces[i].vertex(1)));
                          vert_for_plane_fit.push_back(Point(fixed_faces[i].vertex(0)));
                          original_faces.push_back(Triangle(ref,Point(fixed_faces[i].vertex(1)),Point(fixed_faces[i].vertex(0))));
                   }
              }

          }


        for(int i=0;i<counter;i++)
        {

            test_1 = -1; test_2 = -1,test_3 = -1, test_4 = -1;
            int bf_pointid_1 = -1,bf_pointid_2=-1;

            bool eq_tr = 0;
            for(int t =0;t<original_faces.size();t++)
            {
                if(triEqual(del_face_detail[i].tri,original_faces[t]))
                    eq_tr = 1;
            }
            if(eq_tr)
                continue;

            float queryPoint3[3] = {del_face_detail[i].tri.vertex(0).x(), del_face_detail[i].tri.vertex(0).y(),del_face_detail[i].tri.vertex(0).z()};
            pcdIndex->knnSearch(&queryPoint3[0], vert_indx, &ret_index[0], &out_dist_sqr[0]);
            float queryPoint4[3] = {del_face_detail[i].tri.vertex(1).x(), del_face_detail[i].tri.vertex(1).y(),del_face_detail[i].tri.vertex(1).z()};
            pcdIndex->knnSearch(&queryPoint4[0], vert_indx, &ret_index_1[0], &out_dist_sqr[0]);
            float queryPoint5[3] = {del_face_detail[i].tri.vertex(2).x(), del_face_detail[i].tri.vertex(2).y(),del_face_detail[i].tri.vertex(2).z() };
            pcdIndex->knnSearch(&queryPoint5[0], vert_indx, &ret_index_2[0], &out_dist_sqr[0]);
            buffer.clear();
            if(nbr_tri_of_point[ret_index[0]].done ==13 || nbr_tri_of_point[ret_index_1[0]].done ==13 || nbr_tri_of_point[ret_index_2[0]].done ==13)
            {
                    ;
            }
            else
            {
                    temp_tri.clear();
                    if(ref == del_face_detail[i].tri.vertex(0))
                     {
                          temp_tri.push_back(Triangle(del_face_detail[i].tri.vertex(1),del_face_detail[i].tri.vertex(2),del_face_detail[i].tri.vertex(0)));
                          bf_pointid_1 = ret_index_1[0];
                          temp_tri.push_back(Triangle(del_face_detail[i].tri.vertex(2),del_face_detail[i].tri.vertex(0),del_face_detail[i].tri.vertex(1)));
                          bf_pointid_2 = ret_index_2[0];
                          buffer.clear();
                          if(nbr_tri_buffer[ret_index_1[0]].facelist.size()>0)
                          {
                              for(int k =0;k<nbr_tri_buffer[ret_index_1[0]].facelist.size();k++)
                              {
                                buffer.push_back(nbr_tri_buffer[ret_index_1[0]].facelist[k]);
                              }
                              for(int k =0;k<nbr_tri_of_point[ret_index_1[0]].facelist.size();k++)
                              {
                                buffer.push_back(nbr_tri_of_point[ret_index_1[0]].facelist[k]);
                              }

                              test_1 = test(del_face_detail[i].tri.vertex(1),buffer,del_face_detail[i].tri.vertex(2),del_face_detail[i].tri.vertex(0));
                          }
                          else
                          {
                            buffer = nbr_tri_of_point[ret_index_1[0]].facelist;
                            test_1 = test(del_face_detail[i].tri.vertex(1),buffer,del_face_detail[i].tri.vertex(2),del_face_detail[i].tri.vertex(0));
                          }

                          tri_align_vert.clear();
                          tri_align_vert = buffer;
                          tri_align_vert.push_back(Triangle(del_face_detail[i].tri.vertex(1),del_face_detail[i].tri.vertex(2),del_face_detail[i].tri.vertex(0)));
                          test_3 = check_manifold(del_face_detail[i].tri.vertex(1),tri_align_vert);


                          buffer.clear();
                          if(nbr_tri_buffer[ret_index_2[0]].facelist.size()>0)
                          {

                              for(int k =0;k<nbr_tri_buffer[ret_index_2[0]].facelist.size();k++)
                              {
                                buffer.push_back(nbr_tri_buffer[ret_index_2[0]].facelist[k]);
                              }
                              for(int k =0;k<nbr_tri_of_point[ret_index_2[0]].facelist.size();k++)
                              {
                                buffer.push_back(nbr_tri_of_point[ret_index_2[0]].facelist[k]);
                              }

                              test_2 = test(del_face_detail[i].tri.vertex(2),buffer,del_face_detail[i].tri.vertex(0),del_face_detail[i].tri.vertex(1));
                          }
                          else
                          {
                            buffer = nbr_tri_of_point[ret_index_2[0]].facelist;
                            test_2 = test(del_face_detail[i].tri.vertex(2),buffer,del_face_detail[i].tri.vertex(0),del_face_detail[i].tri.vertex(1));
                          }


                          if(test_1 || test_2)
                              continue;

                          tri_align_vert.clear();
                          tri_align_vert = buffer;
                          tri_align_vert.push_back(Triangle(del_face_detail[i].tri.vertex(2),del_face_detail[i].tri.vertex(0),del_face_detail[i].tri.vertex(1)));
                          test_4 = check_manifold(del_face_detail[i].tri.vertex(2),tri_align_vert);

                          if(test_3 < 0 || test_4 <0)
                              continue;
                          if(test_1 == 0 && test_2 ==0 && test_3 >= 0 && test_4 >=0)
                          {
                              vert_for_plane_fit.push_back(del_face_detail[i].tri.vertex(1));
                              p1 = del_face_detail[i].tri.vertex(1);
                              vert_for_plane_fit.push_back(del_face_detail[i].tri.vertex(2));
			                        p2 = del_face_detail[i].tri.vertex(2);
                          }
                     }
                     else if(ref == del_face_detail[i].tri.vertex(1))
                     {
                        temp_tri.clear();
                        temp_tri.push_back(Triangle(del_face_detail[i].tri.vertex(0),del_face_detail[i].tri.vertex(1),del_face_detail[i].tri.vertex(2)));
                        bf_pointid_1 = ret_index[0];
                        temp_tri.push_back(Triangle(del_face_detail[i].tri.vertex(2),del_face_detail[i].tri.vertex(0),del_face_detail[i].tri.vertex(1)));
                        bf_pointid_2 = ret_index_2[0];
                        buffer.clear();

                       if(nbr_tri_buffer[ret_index[0]].facelist.size()>0)
                       {
                          for(int k =0;k<nbr_tri_buffer[ret_index[0]].facelist.size();k++)
                          {
                            buffer.push_back(nbr_tri_buffer[ret_index[0]].facelist[k]);
                          }
                          for(int k =0;k<nbr_tri_of_point[ret_index[0]].facelist.size();k++)
                          {
                            buffer.push_back(nbr_tri_of_point[ret_index[0]].facelist[k]);
                          }
                         test_1 = test(del_face_detail[i].tri.vertex(0),buffer,del_face_detail[i].tri.vertex(1),del_face_detail[i].tri.vertex(2));
                       }
                       else
                       {
                          buffer = nbr_tri_of_point[ret_index[0]].facelist;
                          test_1 = test(del_face_detail[i].tri.vertex(0),buffer,del_face_detail[i].tri.vertex(1),del_face_detail[i].tri.vertex(2));
                       }


                       tri_align_vert.clear();
                       tri_align_vert = buffer;

                       tri_align_vert.push_back(Triangle(del_face_detail[i].tri.vertex(0),del_face_detail[i].tri.vertex(1),del_face_detail[i].tri.vertex(2)));
                       test_3 = check_manifold(del_face_detail[i].tri.vertex(0),tri_align_vert);

                       buffer.clear();
                       if(nbr_tri_buffer[ret_index_2[0]].facelist.size()>0)
                       {
                          for(int k =0;k<nbr_tri_buffer[ret_index_2[0]].facelist.size();k++)
                          {
                            buffer.push_back(nbr_tri_buffer[ret_index_2[0]].facelist[k]);
                          }
                          for(int k =0;k<nbr_tri_of_point[ret_index_2[0]].facelist.size();k++)
                          {
                            buffer.push_back(nbr_tri_of_point[ret_index_2[0]].facelist[k]);
                          }

                         test_2 = test(del_face_detail[i].tri.vertex(2),buffer,del_face_detail[i].tri.vertex(0),del_face_detail[i].tri.vertex(1));
                       }
                       else
                       {
                         buffer = nbr_tri_of_point[ret_index_2[0]].facelist;
                         test_2 = test(del_face_detail[i].tri.vertex(2),buffer,del_face_detail[i].tri.vertex(0),del_face_detail[i].tri.vertex(1));
                       }


                          if(test_1 || test_2)
                              continue;

                          tri_align_vert.clear();
                          tri_align_vert = buffer;

                          tri_align_vert.push_back(Triangle(del_face_detail[i].tri.vertex(2),del_face_detail[i].tri.vertex(0),del_face_detail[i].tri.vertex(1)));
                          test_4 = check_manifold(del_face_detail[i].tri.vertex(2),tri_align_vert);

                          if(test_3 < 0 || test_4 <0)
                              continue;
                          if(test_1 == 0 && test_2 ==0 && test_3 >= 0 && test_4 >=0)
                          {
                              vert_for_plane_fit.push_back(del_face_detail[i].tri.vertex(0));
			                        p1 = del_face_detail[i].tri.vertex(0);
                              vert_for_plane_fit.push_back(del_face_detail[i].tri.vertex(2));
			                        p2 = del_face_detail[i].tri.vertex(2);
                          }
                     }
                     else
                     {

                       temp_tri.clear();
                       temp_tri.push_back(Triangle(del_face_detail[i].tri.vertex(0),del_face_detail[i].tri.vertex(1),del_face_detail[i].tri.vertex(2)));
                       bf_pointid_1 = ret_index[0];
                       temp_tri.push_back(Triangle(del_face_detail[i].tri.vertex(1),del_face_detail[i].tri.vertex(0),del_face_detail[i].tri.vertex(2)));
                       bf_pointid_2 = ret_index_2[0];
                       buffer.clear();
                       if(nbr_tri_buffer[ret_index[0]].facelist.size()>0)
                       {
                         for(int k =0;k<nbr_tri_buffer[ret_index[0]].facelist.size();k++)
                         {
                           buffer.push_back(nbr_tri_buffer[ret_index[0]].facelist[k]);
                         }
                         for(int k =0;k<nbr_tri_of_point[ret_index[0]].facelist.size();k++)
                         {
                           buffer.push_back(nbr_tri_of_point[ret_index[0]].facelist[k]);
                         }
                         test_1 = test(del_face_detail[i].tri.vertex(0),buffer,del_face_detail[i].tri.vertex(1),del_face_detail[i].tri.vertex(2));
                       }
                       else
                       {
                         buffer = nbr_tri_of_point[ret_index[0]].facelist;
                          test_1 = test(del_face_detail[i].tri.vertex(0),buffer,del_face_detail[i].tri.vertex(1),del_face_detail[i].tri.vertex(2));
                       }

                      tri_align_vert.clear();
                      tri_align_vert = buffer;
                      tri_align_vert.push_back(Triangle(del_face_detail[i].tri.vertex(0),del_face_detail[i].tri.vertex(1),del_face_detail[i].tri.vertex(2)));
                      test_3 = check_manifold(del_face_detail[i].tri.vertex(0),tri_align_vert);

                       if(nbr_tri_buffer[ret_index_1[0]].facelist.size()>0)
                       {
                         buffer.clear();
                         for(int k =0;k<nbr_tri_buffer[ret_index_1[0]].facelist.size();k++)
                         {
                           buffer.push_back(nbr_tri_buffer[ret_index_1[0]].facelist[k]);
                         }
                         for(int k =0;k<nbr_tri_of_point[ret_index_1[0]].facelist.size();k++)
                         {
                           buffer.push_back(nbr_tri_of_point[ret_index_1[0]].facelist[k]);
                         }
                         test_2 = test(del_face_detail[i].tri.vertex(1),buffer,del_face_detail[i].tri.vertex(0),del_face_detail[i].tri.vertex(2));
                       }
                       else
                       {
                         buffer = nbr_tri_of_point[ret_index_1[0]].facelist;
                         test_2 = test(del_face_detail[i].tri.vertex(1),buffer,del_face_detail[i].tri.vertex(0),del_face_detail[i].tri.vertex(2));
                       }


                          if(test_1 || test_2)
                              continue;

                         tri_align_vert.clear();
                         tri_align_vert = buffer;
                          tri_align_vert.push_back(Triangle(del_face_detail[i].tri.vertex(1),del_face_detail[i].tri.vertex(2),del_face_detail[i].tri.vertex(0)));
                          test_4 = check_manifold(del_face_detail[i].tri.vertex(1),tri_align_vert);

                         if(test_3 < 0 || test_4 <0)
                              continue;
                          if(test_1 == 0 && test_2 ==0 && test_3 >= 0 && test_4 >=0)
                          {
                              vert_for_plane_fit.push_back(del_face_detail[i].tri.vertex(0));
                              p1 = del_face_detail[i].tri.vertex(0); bf_pointid_1 = ret_index[0];
                              vert_for_plane_fit.push_back(del_face_detail[i].tri.vertex(1));
			                        p2 = del_face_detail[i].tri.vertex(1);  bf_pointid_2 = ret_index_1[0];

                          }
                      }

                  if(test_1 == 0 && test_2 ==0 && test_3 >= 0 && test_4 >=0)
                  {

		                  linear_least_squares_fitting_3(vert_for_plane_fit.begin(),vert_for_plane_fit.end(),plane,CGAL::Dimension_tag<0>());

                      p_proj1= plane.projection(p1);
                      p_proj2 = plane.projection(p2);
                      p_ref = plane.projection(ref);
                      for(int j=0;j<original_faces.size();j++)
                      {
                            temp_p_proj1= plane.projection(original_faces[j].vertex(1));
                            temp_p_proj2 = plane.projection(original_faces[j].vertex(2));

                            intrsect = check_overlap_trn(Triangle(p_ref,temp_p_proj1,temp_p_proj2),Triangle(p_ref,p_proj1,p_proj2),plane);

                            if(intrsect)
                            {
                                 vert_for_plane_fit.erase(vert_for_plane_fit.end() - 1);
                                 vert_for_plane_fit.erase(vert_for_plane_fit.end() - 1);
                                 break;
                            }
                        }

                        if(!intrsect)
                        {
                            original_faces.push_back(Triangle(ref,p1,p2));
                            additional_faces.push_back(Triangle(ref,p1,p2));
                            nbr_tri_buffer[bf_pointid_1].facelist.push_back(temp_tri[0]);
                            nbr_tri_buffer[bf_pointid_2].facelist.push_back(temp_tri[1]);
                            if(check_manifold(ref,original_faces) == 0)
                               return additional_faces;
                        } // end of if(!res)
                  }
            }
       } // end of for(int i=1;i<counter;i++)*/

     return additional_faces;

  }



//The function is used to get the current time of the day.
//This in tern is used to measure the run time of the program
//get_timestamp ()
static timestamp_t get_timestamp ()
{
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}


/*computes the distance between two points*/
double distance(Point a, Point b)
{
	return sqrt(((a.x() - b.x())*(a.x() - b.x())) + ((a.y() - b.y())*(a.y() - b.y())) + ((a.z() - b.z())*(a.z() - b.z())));
}

double calculateCircumradius(Triangle t)
{

	return distance(CGAL::circumcenter(t), t.vertex(1));
}


/*writes the triangle to the STL file*/
void writeTriangle(Triangle t, std::ofstream &outFile)
{
  outFile << " facet normal 0.0 0.0 0.0" << std::endl;
  outFile << "  outer loop" << std::endl;

  outFile << "   vertex " << t.vertex(0) << std::endl;
  outFile << "   vertex " << t.vertex(1) << std::endl;
  outFile << "   vertex " << t.vertex(2) << std::endl;

  outFile << "  endloop" << std::endl;
  outFile << " endfacet" << std::endl;
}

bool acompare(Del_Faces_Detail lhs, Del_Faces_Detail rhs)
{
    return (lhs.circum_rad < rhs.circum_rad);
}

//This function is not being called . to be removed

std::vector<Triangle> clean_non_manifold_edges(Point ref,std::vector<Point> vert_for_plane_fit,std::vector<Triangle> original_faces)
{
    Plane plane;
    int res =0,cnt=0;
    std::vector<Triangle> indx_1;
    std::vector<int> indx;
    std::vector<Triangle> temp;
    Point c_cent;
    //Cleaning the Hanging Triangles.
    int flag_1 = 13,flag_2 = 13,true_1 = 13;
    while(true_1)
    {
        true_1 =0;
        for(int i=0;i<original_faces.size();i++)
        {
            flag_1 = 13; flag_2 = 13;
            for(int j =0;j<original_faces.size();j++)
            {
                if(i != j)
                {
                    if(original_faces[i].vertex(1) == original_faces[j].vertex(1) || original_faces[i].vertex(1) == original_faces[j].vertex(2))
                    {
                        flag_1 = 0;
                    }
                    if(original_faces[i].vertex(2) == original_faces[j].vertex(1) || original_faces[i].vertex(2) == original_faces[j].vertex(2))
                    {
                        flag_2 = 0;
                    }
                  }
             }
             if(flag_1 == 13 || flag_2 == 13)
             {
                indx_1.push_back(original_faces[i]);

              }
              else
              {
                temp.push_back(original_faces[i]);
              }

          }
          if(temp.size()<original_faces.size())
          {
            true_1=13;
            original_faces.clear();
            original_faces = temp;
            temp.clear();
          }

    }

     if(indx_1.size()>0)
     {
           vert_for_plane_fit.clear();
           for(int j =0;j<original_faces.size();j++)
              {
                  vert_for_plane_fit.push_back(original_faces[j].vertex(0));
                  vert_for_plane_fit.push_back(original_faces[j].vertex(1));
                  vert_for_plane_fit.push_back(original_faces[j].vertex(2));
              }
    }

  //end of Cleaning the Hanging Triangles.

    linear_least_squares_fitting_3(vert_for_plane_fit.begin(),vert_for_plane_fit.end(),plane,CGAL::Dimension_tag<0>());


    for(int i=0;i<original_faces.size();i++)
      for(int j =i+1;j<original_faces.size();j++)
      {
          if(i!=j)
          {
              res = check_overlap_trn(original_faces[i],original_faces[j],plane);
              if(res)
              {
                indx.push_back(i);
                indx.push_back(j);
              }
          }

      }


    for(int k =0;k<indx.size();k++)
      {
        cnt = std::count(indx.begin(),indx.end(), indx[k]);
        if(cnt > 1)
          indx_1.push_back(original_faces[indx[k]]);
      }

    return indx_1;
}

std::vector<Point> hanging_vertices(Point ref,std::vector<Triangle> face_list)
{
  int flag_1 = 13,flag_2 = 13,true_1 = 13;
  std::vector<Point> hang_vert;

      for(int i=0;i<face_list.size();i++)
      {
          flag_1 = 13; flag_2 = 13;
          for(int j =0;j<face_list.size();j++)
          {
              if(i != j)
              {
                  if(face_list[i].vertex(1) == face_list[j].vertex(1) || face_list[i].vertex(1) == face_list[j].vertex(2))
                  {
                      flag_1 = 0;
                  }
                  if(face_list[i].vertex(2) == face_list[j].vertex(1) || face_list[i].vertex(2) == face_list[j].vertex(2))
                  {
                      flag_2 = 0;
                  }
                }
           }
           if(flag_1 == 13)
           {
              hang_vert.push_back(face_list[i].vertex(1));

            }
             if(flag_2 == 13)
            {
              hang_vert.push_back(face_list[i].vertex(2));
            }

        }
    return hang_vert;
}

bool check_one_nbrhd(std::vector<Point> nighbors,Point h_vert)
{
   for(int k =0;k<nighbors.size();k++)
       if(nighbors[k] == h_vert)
           return true;

   return false;
}
// This function is not being called. to be removed
bool check_for_already_processed(std::vector<Triangle> face_indx_rem, std::vector<Triangle> fixed_tri)
{

    for(int k =0;k<face_indx_rem.size();k++)
        for(int j =0;j<fixed_tri.size();j++)
        {
            if(triEqual(face_indx_rem[k],fixed_tri[j]))
               return true;
        }
    return false;
}

void compn(std::vector<Point> p_mfld,std::vector<int> p_mfld_indx,KDTree* pcdIndex,Nbr_Tri_of_Point* nbr_tri_of_point,int indx)
{
  std::vector<Triangle> face_list;
  std::vector<size_t> ret_index(num_nbr+num_nbr);
  std::vector<float> out_dist_sqr(num_nbr+num_nbr);
  int pmfld =0;
  int flgg = 13;
  int point_id_1 = 0;
  std::vector <Point> hang_vert;
  std::vector<Point> nighbors;
  int entered = 0,iteration = 0;;

  while(flgg == 13)
  {
        flgg=0;
        for(int i=0;i<p_mfld.size();i++)
          {
              if(nbr_tri_of_point[p_mfld_indx[i]].done==13)
                 continue;
              pmfld =0;
              float queryPoint11[3] = {p_mfld[i].x(), p_mfld[i].y(),p_mfld[i].z() };
              pcdIndex->knnSearch(&queryPoint11[0], num_nbr, &ret_index[0], &out_dist_sqr[0]);
              point_id_1 = ret_index[0];

              face_list = nbr_tri_of_point[point_id_1].facelist;

              nighbors.clear();
              for(int i=0;i<num_nbr;i++)
                if(ret_index[i] <= indx && ret_index[i] >= 0)
                      nighbors.push_back(vertice[ret_index[i]]);

              pmfld = check_manifold(p_mfld[i] ,face_list);
              if(pmfld == 0)
                  nbr_tri_of_point[point_id_1].done = 13;
              else if(pmfld < 0)
              {
                std::cout<< " TTTTThe Point with negetive manifold is ="<<p_mfld[i] <<" mfld= "<<pmfld<<std::endl; // to be removed

              }
              else //if(pmfld > 0)
              {

                  hang_vert.clear();
                  Triangle triangle;
                  hang_vert = hanging_vertices(p_mfld[i],face_list);
                  if(hang_vert.size() == 2)
                  {
                      triangle = Triangle(p_mfld[i],hang_vert[0],hang_vert[1]);
                      int test_1 =-1,test_2=-1,test_3 =-1;
                      int p_id_2=0,p_id_3=0;
                      std::vector<Triangle> tri_align_vert;
                      face_list.push_back(triangle);
                      test_3 = check_manifold(triangle.vertex(0),face_list);
                      if(test_3<0)
                           continue;
                      float queryPoint22[3] = {triangle.vertex(1).x(), triangle.vertex(1).y(),triangle.vertex(1).z() };
                      pcdIndex->knnSearch(&queryPoint22[0], vert_indx, &ret_index[0], &out_dist_sqr[0]);
                      p_id_2 = ret_index[0];

                      tri_align_vert.clear();
                      tri_align_vert = nbr_tri_of_point[p_id_2].facelist;
                      tri_align_vert.push_back(Triangle(triangle.vertex(1),triangle.vertex(2),triangle.vertex(0)));
                      test_1 = check_manifold(triangle.vertex(1),tri_align_vert);
                      if(test_1<0)
                           continue;

                      float queryPoint23[3] = {triangle.vertex(2).x(), triangle.vertex(2).y(),triangle.vertex(2).z() };
                      pcdIndex->knnSearch(&queryPoint23[0], vert_indx, &ret_index[0], &out_dist_sqr[0]);
                      p_id_3 = ret_index[0];
                      tri_align_vert.clear();
                      tri_align_vert = nbr_tri_of_point[p_id_3].facelist;
                      tri_align_vert.push_back(Triangle(triangle.vertex(2),triangle.vertex(0),triangle.vertex(1)));
                      test_2 = check_manifold(triangle.vertex(2),tri_align_vert);
                      if(test_2<0)
                          continue;

                      if(test_1 ==0 && test_2 ==0 && test_3 ==0)
                      {
                          mesh.push_back(triangle);
                          nbr_tri_of_point[point_id_1].done =13;
                          nbr_tri_of_point[point_id_1].facelist.push_back(triangle);
                          nbr_tri_of_point[p_id_2].facelist.push_back(Triangle(triangle.vertex(1),triangle.vertex(2),triangle.vertex(0)));
                          nbr_tri_of_point[p_id_2].done =13;
                          nbr_tri_of_point[p_id_3].facelist.push_back(Triangle(triangle.vertex(2),triangle.vertex(0),triangle.vertex(1)));
                          nbr_tri_of_point[p_id_3].done =13;
                          flgg = 13;
                      }
                      else if(test_1 >=0 && test_2 >=0 && check_one_nbrhd(nighbors,hang_vert[0]) && check_one_nbrhd(nighbors,hang_vert[1]))
                      {
                          mesh.push_back(triangle);
                          if(test_3 == 0)
                              nbr_tri_of_point[point_id_1].done =13;
                          nbr_tri_of_point[point_id_1].facelist.push_back(triangle);
                          nbr_tri_of_point[p_id_2].facelist.push_back(Triangle(triangle.vertex(1),triangle.vertex(2),triangle.vertex(0)));
                          nbr_tri_of_point[p_id_3].facelist.push_back(Triangle(triangle.vertex(2),triangle.vertex(0),triangle.vertex(1)));
                          flgg = 13;
                        }
                      else //to be removed
                       {
                             entered++;
                       }
                  } // end of hang verttices
              } // end of else

            
              iteration++;
              face_list.clear();
       }//end of  for loop

    }  // end of while

}

int main(int argv, char** argc)
{
    /*Input and output files*/

         if (argv < 3)
         {
             std::cout << "Refer the ReadMe file for furthur details." << std::endl;
             exit(0);
         }
        std::ifstream infile(argc[1]);
        std::ofstream outFile(argc[2]);
        if(argv > 3)
          num_nbr = std::atof(argc[3]);
        else
          num_nbr =16;

        std::cout<< " K  ="<<num_nbr<<std::endl;
        outFile << "solid testsphere" << std::endl;

        PointCloud pcd;
        int indx =0, num_face=0,iteration = 0;
        pcd.pts.resize(20000000);
        int temp=0;
        Point ty; int yt=0;


        std::vector<Point> nbr_vert;
        std::vector<Triangle> faces;
        std::vector<Triangle> computed_faces,temp_original_faces;
        std::vector<Point> p_mfld,n_mfld;
        std::vector<int> p_mfld_indx;


    /*Reading the point cloud from the input file */
        float a, b, c;
        std::cout << "Reading input..."<< std::endl;
        while (infile >> a >> b >> c)
        {
            vertice.push_back(Point(a, b, c));
            pcd.pts[indx].x = a;
            pcd.pts[indx].y = b;
            pcd.pts[indx].z = c;
            indx++;
            if(indx > 19999999)
            {
                std::cout << "The number points is large. Adjust the array size accordingly ..."<< std::endl;
                exit(0);
            }

        }
        infile.close();
        std::cout << "Reading completed..."<< std::endl;

        timestamp_t t0 = get_timestamp();
        Nbr_Tri_of_Point* nbr_tri_of_point = new Nbr_Tri_of_Point[indx];
        Nbr_Tri_of_Point* nbr_tri_buffer = new Nbr_Tri_of_Point[indx];
        std::cout << "Building the KDTree..."<< std::endl;
        std::vector<size_t> ret_index(num_nbr+num_nbr);
        std::vector<float> out_dist_sqr(num_nbr+num_nbr);
        KDTree pcdIndex(3 /*dim*/, pcd, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
        pcdIndex.buildIndex();
        std::cout << "KDTree created "<< std::endl;

        double circum_rad = 0,min_circum_rad = 9999;
        int pointid = 0;

        for (std::vector<Point>::iterator it = vertice.begin() ; it != vertice.end(); ++it)
        {

            Point pt=*it;
            float queryPoint1[3] = {it->x(), it->y(),it->z() };
            pcdIndex.knnSearch(&queryPoint1[0], num_nbr, &ret_index[0], &out_dist_sqr[0]);

            pointid = ret_index[0];

            if(nbr_tri_of_point[pointid].facelist.size()>0)
             {

                 if(check_manifold(pt,nbr_tri_of_point[pointid].facelist) == 0)
                   {
                       nbr_tri_of_point[pointid].done =13;
                       iteration++;
                       continue;
                   }
             }

            for(int i=0;i<num_nbr;i++)
            {
                if(ret_index[i] <= indx && ret_index[i] >= 0)
                {
                    nbr_vert.push_back(vertice[ret_index[i]]);
                    nbr_tri_buffer[ret_index[i]].facelist.clear();
                }
            }

   //computing the delanauy triangulation
            T.clear();
            T.insert(nbr_vert.begin(),nbr_vert.end());
   /*Finding the number of cells in the Delaunay Triangulation */ //-0.275587 -0.154419 0.000004
            int numberOfCells = T.number_of_finite_cells();
            int numberOfEdges = T.number_of_finite_edges();
            int numberOfFaces = T.number_of_finite_facets();

          //  std::cout << "Number of Cells: " << numberOfCells << std::endl; //
          //  std::cout << "Number of edges: " << numberOfEdges << std::endl;//
          //  std::cout << "Number of faces: " << numberOfFaces << std::endl;



         Triangle trii;
         for(Delaunay::Finite_facets_iterator ffi = T.finite_facets_begin(); ffi != T.finite_facets_end(); ffi++)
           {
                 Facet facet = *ffi;
                 trii = T.triangle(facet);
                 if(pt == trii.vertex(0) || pt == trii.vertex(1) || pt == trii.vertex(2))
                 {
                   faces.push_back(trii);
                 }
           }

           int counter =0;
           std::vector<Triangle>::iterator tr_it;
           for(tr_it = faces.begin(); tr_it != faces.end(); tr_it++)
           {
             del_face_detail[counter].tri = *tr_it;
             circum_rad = calculateCircumradius(*tr_it);
             del_face_detail[counter].circum_rad = circum_rad;

             counter++;
           }

           if(counter ==0)
           {
              iteration++;
              continue;
            }
           std::sort(del_face_detail, del_face_detail+counter, acompare); // sort the Facets based on their circum radius

           std::cout<<" iteration = "<<iteration<< " is processing"<<std::endl;
           std::vector<Triangle> face_indx_rem;

           computed_faces = recon(pt,del_face_detail,counter,nbr_tri_of_point[pointid].facelist, &pcdIndex,nbr_tri_of_point,indx,nbr_tri_buffer);

           int mfld = check_manifold(pt,original_faces);


           if(mfld == 0)
           {

              for(int t = 0;t<computed_faces.size();t++)
                  {
                        mesh.push_back(computed_faces[t]);
                        nbr_tri_of_point[pointid].facelist.push_back(computed_faces[t]);

                        float queryPoint2[3] = {computed_faces[t].vertex(1).x(), computed_faces[t].vertex(1).y(),computed_faces[t].vertex(1).z() };
                        pcdIndex.knnSearch(&queryPoint2[0], vert_indx, &ret_index[0], &out_dist_sqr[0]);
                        nbr_tri_of_point[ret_index[0]].facelist.push_back(Triangle(computed_faces[t].vertex(1),computed_faces[t].vertex(2),computed_faces[t].vertex(0)));

                        float queryPoint3[3] = {computed_faces[t].vertex(2).x(), computed_faces[t].vertex(2).y(),computed_faces[t].vertex(2).z() };
                        pcdIndex.knnSearch(&queryPoint3[0], vert_indx, &ret_index[0], &out_dist_sqr[0]);
                        nbr_tri_of_point[ret_index[0]].facelist.push_back(Triangle(computed_faces[t].vertex(2),computed_faces[t].vertex(0),computed_faces[t].vertex(1)));
                  }
                  nbr_tri_of_point[pointid].done =13;
           }
           else if(mfld > 0)
           {
              p_mfld.push_back(pt);
              p_mfld_indx.push_back(pointid);
           }

           std::cout<<" iteration = "<<iteration<< " finished"<<std::endl;

        //Clearing the vectors and vertices vector
            nbr_vert.clear();
            vert_for_plane_fit.clear();
            faces.clear();  //faces1.clear();
            T.clear();
            computed_faces.clear(); //original_faces.clear();

      iteration++;
   } //end of for loop

   compn(p_mfld,p_mfld_indx,&pcdIndex,nbr_tri_of_point,indx);

   timestamp_t t1 = get_timestamp();
   std::cout<<"Size of Mesh = "<<mesh.size()<<std::endl;

   for(int k=0;k<mesh.size();k++)
   {
      writeTriangle(mesh[k], outFile);

   }

   outFile << " endsolid" << std::endl;
   outFile.close();

   std::cout << "Finished Writing STL" << std::endl;


   double secs = (t1 - t0) / 1000000.0L;
   std::cout << "The execution time : " << secs << std::endl;
} // end of main
