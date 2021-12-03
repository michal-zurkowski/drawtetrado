#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#define INF 2e9

// Logical
using Level = int32_t; // 0 - num_levels
using Index = int32_t; // 0/1/2/3
using ID = int32_t;    // Level * 4 + Index;
// Visual
using Position = int32_t; // 0/1/2/3

const std::vector<std::vector<Position>> permutations = {
    {0, 1, 2, 3}, {3, 0, 1, 2}, {2, 3, 0, 1}, {1, 2, 3, 0},
    {2, 1, 0, 3}, {3, 2, 1, 0}, {0, 3, 2, 1}, {1, 0, 3, 2}};

struct Solution {
  std::vector<Position> positions;
  int score = -1;
  std::vector<std::pair<ID, ID>> dangling_edges;
};

Level LevelOf(ID id) { return id / 4; }

// Score is a penalty for the connection. Lower = better.
int Score(const std::pair<Level, Position> &a,
          const std::pair<Level, Position> &b) {
  if (a.first == b.first) {
    // SIMPLE. Simple Connection on the same level.
    return 0;
  } else if (a.second == b.second && std::abs(a.first - b.first) == 1) {
    // SAME_LEVEL. Simple connection up-down 1 level. Best case scenario.
    return 0;
  } else if ((a.second == 2 && b.second == 2) ||
             (a.second == 3 && b.second == 3)) {
    // RIGHT.
    return 3;
  } else if ((a.second == 2 && b.second == 3) ||
             (a.second == 3 && b.second == 2)) {
    // RIGHT_CROSS.
    return 5;
  } else if ((a.second == 0 && b.second == 0) ||
             (a.second == 1 && b.second == 1)) {
    // LEFT.
    return 3;
  } else if ((a.second == 0 && b.second == 1) ||
             (a.second == 1 && b.second == 0)) {
    // LEFT_CROSS.
    return 5;
  } else if ((a.second == 0 && b.second == 3) ||
             (a.second == 3 && b.second == 0)) {
    // FRONT_CROSS.
    return 7;
  } else if ((a.second == 1 && b.second == 2) ||
             (a.second == 2 && b.second == 1)) {
    // BACK_CROSS.
    return 4;
  } else if ((a.second == 1 && b.second == 3) ||
             (a.second == 3 && b.second == 1) ||
             (a.second == 0 && b.second == 2) ||
             (a.second == 2 && b.second == 0)) {
    // FRONT_TO_BACK. // WORST CASE! This connection should not be allowed.
    return 500;
  }

  // Should not reach this part! Make it worst case?
  fprintf(stderr, "Scoring unknown connection type? {%hd %hd} -> {%hd %hd}\n",
          a.first, a.second, b.first, b.second);
  return 500;
}

void UpdateScoreForLevel(const std::vector<ID> &edges, Level level,
                         Solution *solution) {
  auto score = [&](ID start, ID end) {
    return Score({LevelOf(start), solution->positions[start]},
                 {LevelOf(end), solution->positions[end]});
  };
  // Add score for dangling edges ending at given level.
  for (auto &[start, end] : solution->dangling_edges) {
    if (LevelOf(end) == level) {
      solution->score += score(start, end);
    }
  }
  // Remove dangling edges
  auto &d = solution->dangling_edges;
  std::erase_if(d,
                [level](const auto &e) { return LevelOf(e.second) == level; });
  // Add new edges
  for (Index i = 0; i < 4; ++i) {
    ID start = level * 4 + i;
    ID end = edges[start];
    if (LevelOf(end) > level) {
      solution->dangling_edges.emplace_back(start, end);
    } else if (end != -1) {
      solution->score += score(start, end);
    }
  }
}

bool SameLayout(const std::vector<Position> &previous,
                const std::vector<Position> &current,
                const std::vector<Position> &alignments, Level level) {
  // Both should be the same size of 4 elements.
  for (size_t i = 0; i < 4; i++) {
    if (alignments[previous[i] + (level - 1) * 4] !=
        alignments[current[i] + level * 4]) {
      return false;
    }
  }
  return true;
}

Solution Solve(const std::vector<ID> &edges,
               const std::vector<Level> &rotations,
               const std::vector<ID> &alignments) {
  const int num_levels = edges.size() / 4;
  std::vector<Solution> current, next;

  for (const auto &perm : permutations) {
    current.push_back({.positions = perm, .score = 0});
    UpdateScoreForLevel(edges, 0 /* level 0 */, &current.back());
  }

  // current.push_back({.positions = {0, 1, 2, 3}, .score = 0});
  // UpdateScoreForLevel(edges, 0 /* level 0 */, &current.back());
  // current.push_back({.positions = {0, 3, 2, 1}, .score = 0});
  // UpdateScoreForLevel(edges, 0 /* level 0 */, &current.back());

  int best_score = INF;
  for (int level = 1; level < num_levels; ++level) {
    // fprintf(stderr, "Level: %d\n", level);
    // fprintf(stderr, "Solutions: %lu\n\n", current.size());
    next.clear();
    best_score = INF;
    for (const auto &prev_sol : current) {
      for (const auto &perm : permutations) {
        // Is this level in the same rotation group as the previous one?
        if (rotations[level] != -1 &&
            rotations[level] == rotations[level - 1]) {

          if (!SameLayout(
                  {prev_sol.positions.end() - 4, prev_sol.positions.end()},
                  perm, alignments, level)) {
            // Skip this one as it should be rotated the same way as
            // previous level according to tracts.
            continue;
          }
        }

        auto &next_sol = next.emplace_back(prev_sol);
        next_sol.positions.insert(next_sol.positions.end(), perm.begin(),
                                  perm.end());
        UpdateScoreForLevel(edges, level, &next_sol);
        if (best_score > next_sol.score) {
          best_score = next_sol.score;
        }
      }
    }
    // Drop strictly worse solutions
    if (next.size() > 0) {
      int worst_case = best_score + next.back().dangling_edges.size() * 20;
      std::erase_if(
          next, [worst_case](const auto &e) { return e.score > worst_case; });
    }
    std::swap(current, next);
  }

  Solution best_solution;
  best_score = INF;

  for (const auto &solution : current) {
    if (best_score > solution.score) {
      best_score = solution.score;
      best_solution = solution;
    }
  }

  return best_solution;
};

Solution SolveFailsafe(const std::vector<ID> &edges,
                       const std::vector<Level> &rotations,
                       const std::vector<ID> &alignments) {
  Solution result = Solve(edges, rotations, alignments);

  // tracts are not possible to be implemented.
  // Recalculate ignoring tracts.
  if (result.score == -1) {
    fprintf(stderr, "Unable to include tracts, ignoring.\n");
    std::vector<Level> rotations_fixed(edges.size(), -1);
    result = Solve(edges, rotations_fixed, alignments);
  }

  return result;
}

int main() {
  std::vector<ID> edges;
  std::vector<Level> rotations;
  std::vector<ID> alignments;
  int tmp, nucl_num;
  int in_result __attribute__((unused));

  in_result = scanf(" %d ", &nucl_num);
  for (int i = 0; i < nucl_num; ++i) {
    in_result = scanf(" %d ", &tmp);
    edges.push_back(static_cast<ID>(tmp));
  }

  for (int i = 0; i < LevelOf(nucl_num); ++i) {
    in_result = scanf(" %d ", &tmp);
    rotations.push_back(static_cast<Level>(tmp));
  }

  for (int i = 0; i < nucl_num; ++i) {
    in_result = scanf(" %d ", &tmp);
    alignments.push_back(static_cast<ID>(tmp));
  }

  Solution result = Solve(edges, rotations, alignments);

  // tracts are not possible to be implemented.
  // Recalculate ignoring tracts.
  if (result.score == -1) {
    fprintf(stderr, "Unable to include tracts, ignoring.\n");
    for (auto &rotation : rotations) {
      rotation = -1;
    }
    result = Solve(edges, rotations, alignments);
  }
  printf("%hd", result.positions[0]);
  for (size_t i = 1; i < result.positions.size(); ++i) {
    printf(" %hd", result.positions[i]);
  }

  return 0;
}
