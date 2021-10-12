#include <common/pdb.h>
#include <thread>
#include <fstream>
#include <cds/container/fcpriority_queue.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

namespace PMP = CGAL::Polygon_mesh_processing;

struct processed_frame {
  int frame_idx = -1;
  std::vector<Point_3> out_grid_points;

  bool operator<(const processed_frame &other) const
  {
    return frame_idx > other.frame_idx;
  }
};

using priority_queue = cds::container::FCPriorityQueue<processed_frame>;


Surface_Mesh calc_alpha_shape_geometries(const std::vector<Point_3> &points);
Surface_Mesh calc_surface_union_of_balls(const std::vector<Point_3> &points, int radius);

bool g_verbose = true;

auto calc_intersection(const Surface_Mesh &ligand, Surface_Mesh &site) {
  Surface_Mesh intersection;
  auto ligand_cpy = ligand;
  PMP::corefine_and_compute_intersection(ligand_cpy, site, intersection);
  return intersection;
}

processed_frame process_frame(const pdb_frame &f, const CGAL::Surface_mesh<Point_3> &ligand, double site_r)
{
  std::vector<Point_3> points;
  for (const auto &a : f.atoms) {
    points.push_back(a.pos);
  }
  auto site = calc_surface_union_of_balls(points, site_r);

  auto intersection_vol = PMP::volume(calc_intersection(ligand, site));
  auto site_vol = PMP::volume(site);
  auto score = intersection_vol / site_vol;
  std::cout << score << std::endl;

  //

  return {};
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

  parse_pdb(grid_file, queue);

  priority_queue processed_frames;

  std::vector<std::thread> consumers;
  int num_consumers = 1;
  for (int j = 0; j < num_consumers; ++j) {
    consumers.emplace_back([&] {
      bool items_left = false;
      pdb_frame frame;
      do {
        items_left = !queue.producer_done;
        while (queue.frames.try_dequeue(frame)) {
          items_left = true;
          processed_frames.push(process_frame(frame, ligand_mesh, grid_radius));
        }
      } while (items_left ||
               queue.consumers_done.fetch_add(1, std::memory_order_acq_rel) +
                       1 ==
                   num_consumers);
    });
  }

  //producer.join();
  for (auto &c : consumers) {
    c.join();
  }
}