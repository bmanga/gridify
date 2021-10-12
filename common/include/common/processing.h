#ifndef GRIDIFY_PROCESSING_H
#define GRIDIFY_PROCESSING_H

#include <moodycamel/concurrentqueue.h>
#include <cds/container/fcpriority_queue.h>
#include <thread>

#include "common/pdb.h"

using frame_queue = moodycamel::ConcurrentQueue<pdb_frame>;

struct producer_consumer_queue {
  std::atomic_bool producer_done = false;
  std::atomic<int> consumers_done = 0;

  frame_queue frames;
};

template <class T>
struct processed_frame {
  int frame_idx = -1;
  T data;

  bool operator<(const processed_frame &other) const
  {
    return frame_idx > other.frame_idx;
  }
};

template <class Data>
using processed_queue = cds::container::FCPriorityQueue<processed_frame<Data>>;

template <class Data, class Fn>
auto process_frame_loop(int num_consumers, producer_consumer_queue &queue,
                        processed_queue<Data> &processed_frames,
                        Fn &&fn) {
  return [&, num_consumers] {
    bool items_left = false;
    pdb_frame frame;
    do {
      items_left = !queue.producer_done;
      while (queue.frames.try_dequeue(frame)) {
        items_left = true;
        auto data = std::forward<Fn>(fn)(frame);
        processed_frames.push(
            processed_frame{frame.frame_idx, std::move(data)});
      }
    } while (items_left ||
             queue.consumers_done.fetch_add(1, std::memory_order_acq_rel) + 1 ==
                 num_consumers);
  };
}

template <class Data, class Fn>
void process_serialized_results(int num_consumers, producer_consumer_queue &queue,
                                processed_queue<Data> &processed_frames,
                                Fn &&fn)
{
  int top_frame_idx = -1;
  int next_frame = 0;
  while (!processed_frames.empty() ||
         queue.consumers_done != num_consumers + 1) {
    if (processed_frames.empty()) {
      std::this_thread::sleep_for(std::chrono::milliseconds(50));
      continue;
    }
    do {
      processed_frames.apply([&](const auto &queue) {
        if (!queue.empty()) {
          top_frame_idx = queue.top().frame_idx;
        }
      });
    } while (top_frame_idx != next_frame &&
             (queue.consumers_done != num_consumers + 1));
    ++next_frame;
    processed_frame<Data> pf;

    if (processed_frames.pop(pf)) {
      std::forward<Fn>(fn)(pf.frame_idx, pf.data);
    }
  }
}

#endif  // GRIDIFY_PROCESSING_H
