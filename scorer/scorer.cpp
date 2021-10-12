#include <common/pdb.h>
#include "common/processing.h"
#include <thread>
#include <fstream>
#include <cds/container/fcpriority_queue.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <fmt/format.h>

namespace PMP = CGAL::Polygon_mesh_processing;

struct processed_data {
  double intersection_over_union;
  double intersection_over_ligand;
  double intersection_over_site;
};

Surface_Mesh calc_alpha_shape_geometries(const std::vector<Point_3> &points);
Surface_Mesh calc_surface_union_of_balls(const std::vector<Point_3> &points, int radius);

bool g_verbose = true;

auto calc_intersection(const Surface_Mesh &ligand, Surface_Mesh &site) {
  Surface_Mesh intersection;
  auto ligand_cpy = ligand;
  PMP::corefine_and_compute_intersection(ligand_cpy, site, intersection);

  return intersection;
}

processed_data process_frame(const pdb_frame &f, const CGAL::Surface_mesh<Point_3> &ligand, double site_r)
{
  std::vector<Point_3> points;
  for (const auto &a : f.atoms) {
    points.push_back(a.pos);
  }
  auto site = calc_surface_union_of_balls(points, site_r);

  {
    std::ofstream f ("site.off");
    f << site;
    f.close ();
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
      .intersection_over_union = intersection_vol / union_vol,
      .intersection_over_ligand = intersection_vol / ligand_vol,
      .intersection_over_site = intersection_vol / site_vol
  };
}

int main(int argc, char *argv[])
{
  auto ligand_file = std::ifstream(argv[1]);

  auto grid_file = std::ifstream(argv[2]);

  double grid_radius = parse_pdb_gridify_spacing(grid_file);

  grid_file.seekg(0);

  producer_consumer_queue queue_ligand, queue;


  //std::thread producer([&] { parse_pdb(pdb_file, queue); });

  parse_pdb(ligand_file, queue_ligand);

  pdb_frame frame;
  if (!queue_ligand.frames.try_dequeue(frame)) {
    // TODO error.
  }

  std::vector<Point_3> points;
  for (const auto &a : frame.atoms) {
    points.push_back(a.pos);
  }
  auto ligand_mesh = calc_alpha_shape_geometries(points);

  {
    std::ofstream f("ligand.off");
    f << ligand_mesh;
    f.close();
  }

  parse_pdb(grid_file, queue);

  processed_queue<processed_data> processed_frames;

  std::vector<std::thread> consumers;
  int num_consumers = 1;
  for (int j = 0; j < num_consumers; ++j) {
    consumers.emplace_back(process_frame_loop(num_consumers, queue, processed_frames,
                                              [&](const pdb_frame &frame) {
                                                return process_frame(frame, ligand_mesh, grid_radius);
                                              }));
  }

  process_serialized_results(num_consumers, queue, processed_frames,
                             [&](int frame_idx, processed_data &data) {
                               fmt::print("-- writing frame {}\n", frame_idx);
                               fmt::print("IOU: {:.3}\n", data.intersection_over_union);
                               fmt::print("IOL: {:.3}\n", data.intersection_over_ligand);
                               fmt::print("IOS: {:.3}\n", data.intersection_over_site);
                             });

  //producer.join();
  for (auto &c : consumers) {
    c.join();
  }
}