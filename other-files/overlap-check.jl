
# MC-CF checking if too close or not moved
for k in 1:num_parts
				#println(index,", ",i)
				for q in 1:qpart[1]
					escape_attempts = 0
					while abs(running_config[k] - qpart[2][q]) <= lstar*0.001
						if escape_attempts == 0
							println("Particle $k too close to QPart $q")
						end
						get_out = acc_rej_move(running_config,n,p,k,step_size,qpart,log_form)
						running_config = get_out[1]
						wavefunc = get_out[3]
						escape_attempts += 1
					end
					if escape_attempts > 0
						println("Particle $k Away after $escape_attempts attempts")
					end
				end
				#
				if i > 100 && abs(time_config[k,index-101+samp_freq] - running_config[k]) <= lstar*0.001
					got_away = 0
					escape_attempts = 0
					println("Particle $k Stuck")
					while got_away < 0.5
						get_out = acc_rej_move(running_config,n,p,k,step_size,qpart,log_form)
						running_config = get_out[1]
						wavefunc = get_out[3]
						escape_attempts += 1
						got_away = get_out[2]
					end
					println("Particle $k Free after $escape_attempts attempts")
				end
			end
