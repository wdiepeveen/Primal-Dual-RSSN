export dbg, dbglevel, dbglevel_atleast

using CPUTime
using Printf
using Logging

function dbg(level::Int, message::String)

  global debug_level
  global debug_timer_start
  global io

	if !@isdefined(debug_level) || (level <= debug_level)
		if @isdefined(debug_timer_start)
			if level == 0
				@printf io "%s" message * "\n"
				@printf "%s" message * "\n"
			else
				@printf io "[%09.4f] %s" (CPUtime_us()-debug_timer_start)/1e6 message * "\n"
				@printf "[%09.4f] %s" (CPUtime_us()-debug_timer_start)/1e6 message * "\n"
			end
		end
	end

end

function dbglevel(level::Int)
  global debug_level
  debug_level = level
end

function dbglevel_atleast(level::Int)
	global debug_level
	return (!@isdefined(debug_level) || (level <= debug_level))
end
