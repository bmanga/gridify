#include <core/pdb.h>
#include "core/processing.h"
#include "core/radius.h"
#include "core/pca.h"
#include <thread>
#include <fstream>
#include <cds/container/fcpriority_queue.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <fmt/format.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Polygon_mesh_processing/transform.h>

#define ALL_SETTINGS            \
  X(bool, pca_align, false, "align the grid using PCA, with maximum variance on the X axis") \
  X(std::string, ligand_file, "", "The file of the ligand pdb")                                   \
  X(std::string, grid_file, "", "A file generated by gridify")                               \
  X(double, scale_radius, 1.0, "The scale factor for atomic radii (default 1)")

#include "core/cmdline.inc"

namespace PMP = CGAL::Polygon_mesh_processing;

struct processed_data {
  double site_volume;
  double ligand_volume;
  double intersection_volume;
  double union_volume;
};

Surface_Mesh calc_alpha_shape_geometries(const std::vector<Point_3> &points);
Surface_Mesh calc_surface_union_of_balls(const std::vector<Point_3> &points, double radius);
Surface_Mesh calc_surface_union_of_balls(const std::vector<Point_3> &points, const std::vector<double> &radii);

bool g_verbose = true;

auto calc_intersection(const Surface_Mesh &ligand, Surface_Mesh &site) {
  Surface_Mesh intersection;
  auto ligand_cpy = ligand;
  PMP::corefine_and_compute_intersection(ligand_cpy, site, intersection);

  return intersection;
}

processed_data process_frame(const config &c, const pdb_frame &f, const CGAL::Surface_mesh<Point_3> &ligand, double site_r)
{
  std::vector<Point_3> points;
  for (const auto &a : f.atoms) {
    points.push_back(a.pos);
  }

  if (c.pca_align) {
    points = pca_aligned_points(points);
  }

  auto site = calc_surface_union_of_balls(points, site_r);

  {
    std::ofstream f ("site.off");
    f << site;
    f.close ();
  }

  if (PMP::does_self_intersect(site)) {
    fmt::print("site mesh contains self intersections");
    std::exit(-1);
  }

  auto intersection = calc_intersection(ligand, site);
  {
    std::ofstream f("intersection.off");
    f << intersection;
    f.close();
  }
  auto intersection_vol = PMP::volume(intersection);
  auto site_vol = PMP::volume(site);
  auto ligand_vol = PMP::volume(ligand);
  auto union_vol = ligand_vol + site_vol - intersection_vol;

  return {
    .site_volume = site_vol,
    .ligand_volume = ligand_vol,
    .intersection_volume = intersection_vol,
    .union_volume = union_vol
  };
}

auto gen_ligand_geometry(const config &c, const pdb_frame &f)
{
  std::vector<Point_3> points;
  std::vector<double> radii;
  const auto &radmatch = radius_matcher::get();
  for (const auto &a : f.atoms) {
    points.push_back(a.pos);
    radii.push_back(radmatch.radius(a) * c.scale_radius);
  }

  if (c.pca_align) {
    points = pca_aligned_points(points);
  }
  return calc_surface_union_of_balls(points, radii);
}

int main(int argc, char *argv[])
{
  auto config = parse_cmd_line_settings("scorer",
                                        "calculate intersection between a ligand and a grid file",
                                        argc, argv, [](cxxopts::Options &){});

  auto ligand_file = std::ifstream(config.ligand_file);
  auto grid_file = std::ifstream(config.grid_file);

  double grid_radius = parse_pdb_gridify_spacing(grid_file);

  grid_file.seekg(0);

  producer_consumer_queue queue_ligand, queue;


  //std::thread producer([&] { parse_pdb(pdb_file, queue); });

  parse_pdb(ligand_file, queue_ligand);

  pdb_frame frame;
  if (!queue_ligand.frames.try_dequeue(frame)) {
    // TODO error.
  }

  auto ligand_mesh = gen_ligand_geometry(config, frame);

  if (PMP::does_self_intersect(ligand_mesh)) {
    fmt::print("ligand mesh contains self intersections");
    std::exit(-1);
  }

  {
    std::ofstream f("ligand.off");
    f << ligand_mesh;
    f.close();
  }

  std::thread producer([&] {parse_pdb(grid_file, queue);});

  processed_queue<processed_data> processed_frames;

  std::vector<std::thread> consumers;
  int num_consumers = 1;
  for (int j = 0; j < num_consumers; ++j) {
    consumers.emplace_back(process_frame_loop(num_consumers, queue, processed_frames,
                                              [&](const pdb_frame &frame) {
                                                return process_frame(config, frame, ligand_mesh, grid_radius);
                                              }));
  }


  process_serialized_results(
      num_consumers, queue, processed_frames,
      [&](int frame_idx, processed_data &data) {
        fmt::print("-- writing frame {}\n", frame_idx);
        fmt::print("SITE VOLUME:         {:.3}\n", data.site_volume);
        fmt::print("LIGAND VOLUME:       {:.3}\n", data.ligand_volume);
        fmt::print("INTERSECTION VOLUME: {:.3}\n", data.intersection_volume);
        fmt::print("UNION VOLUME:        {:.3}\n", data.union_volume);
        fmt::print("IOU: {:.3}\n", data.intersection_volume / data.union_volume);
        fmt::print("IOL: {:.3}\n", data.intersection_volume / data.ligand_volume);
        fmt::print("IOS: {:.3}\n", data.intersection_volume / data.site_volume);
      });


  producer.join();
  for (auto &c : consumers) {
    c.join();
  }
}