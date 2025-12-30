// -*- lsst-c++ -*-
/*
 * Zaki's Common Library
 * See License file at the top of the source tree.
 *
 * Copyright (c) 2023 Mohammadreza Zakeri
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef Zaki_Util_Instrumentor_H
#define Zaki_Util_Instrumentor_H

/**
 * @file Instrumentor.hpp
 * @brief Lightweight Chrome Tracing (chrome://tracing) profiler.
 *
 * Design goals:
 *  - Thread-safe session lifecycle (Begin/End vs writes).
 *  - Low overhead per event (avoid disk I/O under contention).
 *  - Reliable timestamps via std::chrono::steady_clock with session-relative time base.
 *
 * Output format:
 *  - JSON traceEvents compatible with Chrome Tracing / Perfetto UI.
 *
 * Usage:
 *  @code
 *  Zaki::Util::Instrumentor::BeginSession("MyRun", "results.json");
 *  {
 *      PROFILE_SCOPE("HotScope");
 *      // ...
 *  }
 *  Zaki::Util::Instrumentor::EndSession();
 *  @endcode
 */

#include <atomic>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <mutex>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include "Zaki/String/Directory.hpp"

namespace Zaki::Util
{

/**
 * @brief Singleton instrumentation backend that writes Chrome trace JSON.
 *
 * Thread-safety model:
 *  - Session lifecycle and file stream are protected by an internal mutex.
 *  - Profiling events are accumulated per-thread in a thread_local buffer.
 *  - Buffers are flushed to the shared stream under the same mutex.
 *
 * Important note:
 *  - EndSession() flushes the *calling thread's* buffer and closes the stream.
 *    For complete reliability, call EndSession() after worker threads have stopped
 *    emitting profiling events (e.g., after joining threads), or call Flush() on
 *    each worker before it exits.
 */
class Instrumentor
{
	friend class InstrumentationTimer; // <-- allow timer fast-path access
  public:
	/** @brief Start a profiling session. If a session is already active, it is ended first. */
	static void BeginSession(const std::string &sessionName,
							 const Zaki::String::Directory &filepath = "results.json");

	/** @brief End the active profiling session (no-op if none). */
	static void EndSession();

	/**
	 * @brief Flush the calling thread's buffered events to the output stream.
	 *
	 * Use this in worker threads if they finish before the main thread calls EndSession().
	 */
	static void Flush();

	/**
	 * @brief Low-level event record (called by timers/macros).
	 *
	 * @param name      Event name (scope/function). Only needs to remain valid for the call duration.
	 * @param start_us  Start timestamp in microseconds relative to session start.
	 * @param dur_us    Duration in microseconds.
	 * @param tid       Thread identifier (hashed).
	 */
	static void RecordEvent(std::string_view name,
							std::int64_t start_us,
							std::int64_t dur_us,
							std::uint64_t tid);

	/** @brief Retrieve singleton instance. */
	static Instrumentor &Get();

	/** @brief Destructor ends an active session (best-effort). */
	~Instrumentor();

	Instrumentor(const Instrumentor &) = delete;
	Instrumentor &operator=(const Instrumentor &) = delete;

  private:
	Instrumentor() = default;

	// Internal session control (must be called via static wrappers).
	void IBeginSession(const std::string &sessionName,
					   const Zaki::String::Directory &filepath);
	void IEndSession();

	// Flush calling thread's buffer (internal; static Flush() calls this).
	void IFlushThreadBuffer();

	// JSON helpers (require mutex held because they touch the stream and counters).
	void WriteHeaderLocked();
	void WriteFooterLocked();
	void WriteEventsLocked(const std::vector<std::string> &events);

	// Build a single JSON event string (no locking; purely local work).
	static std::string BuildEventJSON(std::string_view name,
									  std::int64_t start_us,
									  std::int64_t dur_us,
									  std::uint64_t tid);

	// Session state (protected by mutex for stream + counters; active flag is atomic for fast path).
	std::string m_sessionName{"None"};
	std::ofstream m_output;
	std::mutex m_mutex;

	std::atomic<bool> m_active{false};

	bool IsActiveFast() const noexcept
	{
		return m_active.load(std::memory_order_acquire);
	}

	// Session timebase.
	std::chrono::steady_clock::time_point m_sessionStart{};

	// Comma management: number of events already written to file (global).
	std::uint64_t m_eventCount{0};

	// Larger ofstream buffer to reduce syscalls.
	std::vector<char> m_streamBuffer;

	// Per-thread buffered events (thread_local).
	struct ThreadBuffer
	{
		std::vector<std::string> events;
		std::size_t approxBytes{0};

		~ThreadBuffer(); // flush on thread exit

		void clear()
		{
			events.clear();
			approxBytes = 0;
		}
	};

	// Accessor for thread_local buffer for this Instrumentor instance.
	static ThreadBuffer &TLSBuffer();

	// Flush threshold (per thread). Tunable.
	static constexpr std::size_t kFlushThresholdBytes = 256 * 1024; // 256 KB
};

/**
 * @brief RAII timer that records a profiling event on destruction.
 *
 * Prefer constructing it via PROFILE_SCOPE / PROFILE_FUNCTION macros.
 */
class InstrumentationTimer
{
  public:
	/** @brief Start timing a named scope. */
	explicit InstrumentationTimer(std::string_view name);

	/** @brief If not stopped, records the event. */
	~InstrumentationTimer();

	InstrumentationTimer(const InstrumentationTimer &) = delete;
	InstrumentationTimer &operator=(const InstrumentationTimer &) = delete;

	/** @brief Stop and record event (idempotent). */
	void Stop();

  private:
	std::string_view m_name;
	std::chrono::steady_clock::time_point m_start;
	bool m_stopped{false};
};

/**
 * @brief RAII session helper (optional convenience).
 *
 * @code
 *  Zaki::Util::ScopedProfileSession s("Run", "results.json");
 *  // ...
 * @endcode
 */
class ScopedProfileSession
{
  public:
	ScopedProfileSession(const std::string &sessionName,
						 const Zaki::String::Directory &filepath = "results.json")
	{
		Instrumentor::BeginSession(sessionName, filepath);
	}
	~ScopedProfileSession()
	{
		Instrumentor::EndSession();
	}

	ScopedProfileSession(const ScopedProfileSession &) = delete;
	ScopedProfileSession &operator=(const ScopedProfileSession &) = delete;
};

} // namespace Zaki::Util

//==============================================================
// Profiling Macros (enable/disable cleanly at build time)
//==============================================================

// Token pasting helpers (no external __CONCAT dependency).
#define ZAKI_DETAIL_CONCAT_INNER(a, b) a##b
#define ZAKI_DETAIL_CONCAT(a, b) ZAKI_DETAIL_CONCAT_INNER(a, b)

#if defined(_MSC_VER)
#define Z_PRETTY_FUNCTION __FUNCSIG__
#else
#define Z_PRETTY_FUNCTION __PRETTY_FUNCTION__
#endif

/**
 * Define ZAKI_ENABLE_PROFILING=1 at compile time to enable.
 * Example (CMake):
 *   target_compile_definitions(MyTarget PRIVATE ZAKI_ENABLE_PROFILING=1)
 */
#ifndef ZAKI_ENABLE_PROFILING
#define ZAKI_ENABLE_PROFILING 1
#endif

#if ZAKI_ENABLE_PROFILING
#define PROFILE_SCOPE(name) \
	::Zaki::Util::InstrumentationTimer ZAKI_DETAIL_CONCAT(_zaki_timer_, __COUNTER__)(name)

#define PROFILE_FUNCTION() PROFILE_SCOPE(Z_PRETTY_FUNCTION)
#else
#define PROFILE_SCOPE(name) ((void)0)
#define PROFILE_FUNCTION() ((void)0)
#endif

#endif /*Zaki_Util_Instrumentor_H*/