// Copyright (c) 2013-2018 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
// Copyright (c) 2013-2018 Centre national de la recherche scientifique (CNRS)
// Copyright (c) 2018-2023 Simons Foundation
// Copyright (c) 2017 Hugo U.R. Strand
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You may obtain a copy of the License at
//     https://www.gnu.org/licenses/gpl-3.0.txt
//
// Authors: Michel Ferrero, Henri Menke, Olivier Parcollet, Priyanka Seth, Hugo U. R. Strand, Nils Wentzell, Thomas Ayral

#include "./concepts.hpp"
#include "./mc_generic.hpp"
#include "../utility/signal_handler.hpp"
#include "../utility/timestamp.hpp"

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <mpi/mpi.hpp>
#include <mpi/monitor.hpp>
#include <mpi/vector.hpp>
#include <nda/macros.hpp>

#include <algorithm>
#include <complex>
#include <cstdint>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace triqs::mc_tools {

  template <DoubleOrComplex MCSignType> int mc_generic<MCSignType>::run(run_param_t const &params) {
    EXPECTS(params.cycle_length > 0);
    EXPECTS(params.stop_callback);
    EXPECTS(params.after_cycle_duty);

    // set the initial sign if specified
    if (params.initial_sign != default_initial_sign) sign_ = params.initial_sign;

    // reset move statistics
    moves_.clear_statistics();

    // start timer
    run_timer_ = {};
    run_timer_.start();

    // start signal handler
    triqs::signal_handler::start();

    // start MPI monitors
    std::unique_ptr<mpi::monitor> exception_monitor;
    if (params.propagate_exception and mpi::has_env) exception_monitor = std::make_unique<mpi::monitor>(params.comm);
    std::unique_ptr<mpi::monitor> cycle_monitor;
    if (params.continue_after_ncycles_done and mpi::has_env) cycle_monitor = std::make_unique<mpi::monitor>(params.comm);

    // prepare the simulation parameters and statistics
    auto const rank            = params.comm.rank();
    bool stop_flag             = params.ncycles == 0;
    std::int64_t cycle_counter = 1;
    double next_print_info     = 0.1;
    double next_check_except   = params.check_exception_interval;
    double next_check_cycles   = params.check_cycles_interval;
    percentage_done_           = stop_flag ? 100 : 0;
    nmeasures_done_            = 0;

    // run simulation
    for (; !stop_flag; ++cycle_counter) {
      try {
        // do cycle_length MC steps / cycle
        for (std::int64_t i = 0; i < params.cycle_length; i++) {
          if (triqs::signal_handler::received()) throw triqs::signal_handler::exception{};
          // Metropolis step
          double r = moves_.attempt();
          if (rng_() < std::min(1.0, r)) {
            sign_ *= moves_.accept();
          } else {
            moves_.reject();
          }
          ++config_id_;
        }
        // after cycle duties
        params.after_cycle_duty();
        if (params.enable_calibration) moves_.calibrate(params.comm);
        if (params.enable_measures) {
          nmeasures_done_++;
          for (auto &m : measures_aux_) m();
          measures_.accumulate(sign_);
        }
      } catch (triqs::signal_handler::exception const &) {
        // current cycle is interrupted, simulation is stopped below
        std::cerr << fmt::format("[Rank {}] Signal caught in mc_generic::run: Stopping the simulation.\n", rank);
      } catch (std::exception const &err) {
        // log the exception and node number, either abort or report to other nodes
        std::cerr << fmt::format("[Rank {}] Error int mc_generic::run: Exception occured:\n{}\n", rank, err.what());
        if (params.propagate_exception) {
          exception_monitor->report_local_event();
          break;
        } else {
          params.comm.abort(2);
        }
      }

      // recompute fraction done and runtime so far
      percentage_done_ = cycle_counter * 100.0 / params.ncycles;
      double runtime   = run_timer_;

      // is it time to print simulation info
      if (runtime > next_print_info) {
        // increase time interval non-linearly until next print and print info
        next_print_info = 1.25 * runtime + 2.0;
        print_sim_info(params, cycle_counter);
      }

      // is it time to check for exceptions on other ranks
      if (exception_monitor && runtime > next_check_except) {
        // increase time until next check and check other ranks
        next_check_except += params.check_exception_interval;
        stop_flag |= exception_monitor->event_on_any_rank();
      }

      // have we done all requested cycles
      if (percentage_done_ >= 100) {
        if (cycle_monitor) {
          // if continue_after_ncycles_done == true, report to other ranks
          cycle_monitor->report_local_event();
          // is it time to check if the other ranks are finished as well
          if (runtime > next_check_cycles) {
            // increase time until next check and check other ranks
            next_check_cycles += params.check_cycles_interval;
            stop_flag |= cycle_monitor->event_on_all_ranks();
          }
        } else {
          // if continue_after_ncycles_done == false, stop the simulation
          stop_flag = true;
        }
      }

      // update stop flag
      stop_flag |= (params.stop_callback() || triqs::signal_handler::received());
    }

    // stop timer
    run_timer_.stop();

    // update statistics
    --cycle_counter;
    ncycles_done_ += cycle_counter;

    // print final simulation info
    print_sim_info(params, cycle_counter);

    // stop signal handler
    int status = (percentage_done_ >= 100 ? 0 : (triqs::signal_handler::received() ? 2 : 1));
    triqs::signal_handler::stop();

    // stop MPI monitors
    if (exception_monitor) {
      exception_monitor->finalize_communications();
      if (exception_monitor->event_on_any_rank())
        throw std::runtime_error(fmt::format("[Rank {}] MC simulation stopped because an exception occurred on one of the MPI ranks", rank));
    }
    if (cycle_monitor) cycle_monitor->finalize_communications();

    // final reports
    if (status == 1) report_ << fmt::format("[Rank {}] MC simulation stopped because stop_callback() returned true\n", rank);
    if (status == 2) report_ << fmt::format("[Rank {}] MC simulation stopped because a signal has been received\n", rank);

    return status;
  }

  template <DoubleOrComplex MCSignType> int mc_generic<MCSignType>::warmup(run_param_t const &params) {
    report_(3) << fmt::format("[Rank {}] Performing warum up phase...\n", params.comm.rank());
    auto p            = params;
    p.enable_measures = false;
    auto status       = run(p);
    warmup_timer_     = run_timer_;
    return status;
  }

  template <DoubleOrComplex MCSignType> int mc_generic<MCSignType>::accumulate(run_param_t const &params) {
    report_(3) << fmt::format("[Rank {}] Performing accumulation phase...\n", params.comm.rank());
    auto p               = params;
    p.enable_measures    = true;
    p.enable_calibration = false;
    auto status          = run(p);
    acc_timer_           = run_timer_;
    return status;
  }

  template <DoubleOrComplex MCSignType>
  int mc_generic<MCSignType>::run(std::int64_t ncycles, std::int64_t cycle_length, std::function<bool()> stop_callback, bool enable_measures,
                                  mpi::communicator c, bool enable_calibration) {
    return run({.ncycles            = ncycles,
                .cycle_length       = cycle_length,
                .stop_callback      = stop_callback,
                .comm               = c,
                .enable_measures    = enable_measures,
                .enable_calibration = enable_calibration});
  }

  template <DoubleOrComplex MCSignType>
  int mc_generic<MCSignType>::warmup(std::int64_t ncycles, std::int64_t cycle_length, std::function<bool()> stop_callback, MCSignType initial_sign,
                                     mpi::communicator c) {
    return warmup({.ncycles = ncycles, .cycle_length = cycle_length, .stop_callback = stop_callback, .initial_sign = initial_sign, .comm = c});
  }

  template <DoubleOrComplex MCSignType>
  int mc_generic<MCSignType>::warmup(std::int64_t ncycles, std::int64_t cycle_length, std::function<bool()> stop_callback, mpi::communicator c) {
    return warmup({.ncycles = ncycles, .cycle_length = cycle_length, .stop_callback = stop_callback, .comm = c});
  }

  template <DoubleOrComplex MCSignType>
  int mc_generic<MCSignType>::accumulate(std::int64_t ncycles, std::int64_t cycle_length, std::function<bool()> stop_callback, mpi::communicator c) {
    return accumulate({.ncycles = ncycles, .cycle_length = cycle_length, .stop_callback = stop_callback, .comm = c});
  }

  template <DoubleOrComplex MCSignType>
  int mc_generic<MCSignType>::warmup_and_accumulate(std::int64_t ncycles_warmup, std::int64_t ncycles_acc, std::int64_t cycle_length,
                                                    std::function<bool()> stop_callback, MCSignType initial_sign, mpi::communicator c) {
    auto status =
       warmup({.ncycles = ncycles_warmup, .cycle_length = cycle_length, .stop_callback = stop_callback, .initial_sign = initial_sign, .comm = c});
    if (status == 0) status = accumulate({.ncycles = ncycles_acc, .cycle_length = cycle_length, .stop_callback = stop_callback, .comm = c});
    return status;
  }

  template <DoubleOrComplex MCSignType>
  int mc_generic<MCSignType>::warmup_and_accumulate(std::int64_t ncycles_warmup, std::int64_t ncycles_acc, std::int64_t cycle_length,
                                                    std::function<bool()> stop_callback, mpi::communicator c) {
    return warmup_and_accumulate(ncycles_warmup, ncycles_acc, cycle_length, stop_callback, default_initial_sign, c);
  }

  template <DoubleOrComplex MCSignType> void mc_generic<MCSignType>::collect_results(mpi::communicator const &c) {
    report_(3) << fmt::format("[Rank {}] Collect results: Waiting for all MPI processes to finish accumulating...\n", c.rank());

    // collect results from all MPI processes
    measures_.collect_results(c);
    moves_.collect_statistics(c);
    auto tot_nmeasures = mpi::reduce(nmeasures_done_, c);
    auto tot_duration  = mpi::reduce(get_accumulation_time(), c);

    // generate rank dependent output string
    std::string info{"\n"};
    info += fmt::format("[Rank {}] Warmup duration: {:.4f} seconds [{}]\n", c.rank(), get_warmup_time(), get_warmup_time_HHMMSS());
    info += fmt::format("[Rank {}] Simulation duration: {:.4f} seconds [{}]\n", c.rank(), get_accumulation_time(), get_accumulation_time_HHMMSS());
    info += fmt::format("[Rank {}] Number of measures: {}\n", c.rank(), nmeasures_done_);
    info += fmt::format("[Rank {}] Cycles (measures) / second: {:.2e}\n", c.rank(), nmeasures_done_ / get_accumulation_time());
    info += fmt::format("[Rank {}] Measurement durations:\n{}", c.rank(), measures_.get_timings(fmt::format("[Rank {}]   ", c.rank())));
    info += fmt::format("[Rank {}] Move statistics:\n{}", c.rank(), moves_.get_statistics(fmt::format("[Rank {}]   ", c.rank())));

    // gather all output strings on rank 0 to print in order
    auto all_infos_vec = mpi::gather(std::vector<char>{info.begin(), info.end()}, c);
    if (c.rank() == 0) {
      std::string all_infos{all_infos_vec.begin(), all_infos_vec.end()};
      report_(3) << all_infos;
      std::string more_info{"\n"};
      more_info += fmt::format("Total number of measures: {}\n", tot_nmeasures);
      more_info += fmt::format("Total cycles (measures) / second: {:.2e}\n", tot_nmeasures / tot_duration);
      report_(2) << more_info;
    }
  }

  template <DoubleOrComplex MCSignType> void mc_generic<MCSignType>::print_sim_info(run_param_t const &params, std::int64_t cycle_counter) {
    // current simulation parameters
    auto const rank           = params.comm.rank();
    double const runtime      = run_timer_;
    auto const cycles_per_sec = cycle_counter / runtime;

    // do the printing
    if (percentage_done_ < 0) {
      report_(3) << fmt::format("[Rank {}] {} cycle {}, {:.2e} cycles/sec\n", rank, utility::timestamp(), cycle_counter, cycles_per_sec);
    } else {
      report_(3) << fmt::format("[Rank {}] {} {:>6.2f}% done, ETA {}, cycle {} of {}, {:.2e} cycles/sec\n", rank, utility::timestamp(),
                                percentage_done_, utility::estimate_time_left(params.ncycles, cycle_counter, run_timer_), cycle_counter,
                                params.ncycles, cycles_per_sec);
    }
    if (params.enable_measures) report_(3) << measures_.report();
  }

  // Explicit template instantiations.
  template class mc_generic<double>;
  template class mc_generic<std::complex<double>>;

} // namespace triqs::mc_tools
