
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2017 by Hugo U.R. Strand
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <chrono>

namespace triqs {
  namespace utility {

    class timer {
      public:
      using clock_t = std::chrono::high_resolution_clock;

      private:
      clock_t::time_point start_time;
      clock_t::duration total_time = clock_t::duration(0);
      bool running = false;

      public:
      void start() {
        running    = true;
        start_time = clock_t::now();
      }
      void stop() {
        total_time += clock_t::now() - start_time;
        running = false;
      }
      bool is_running() const { return running; }
      operator double() const {
        std::chrono::duration<double> total_time_seconds(total_time);
        if (is_running()) total_time_seconds += clock_t::now() - start_time;
        return total_time_seconds.count();
      }
    };
  }
}
