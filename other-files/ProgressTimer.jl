

if !@isdefined(TabLevel)
    TabLevel = ""
end
println(TabLevel*"Open ProgressTimer.jl")
TabLevel=TabLevel*"    "

using Dates


const sec_len=Dates.Millisecond(Dates.Second(1))
const min_len=Dates.Millisecond(Dates.Minute(1))
const hr_len=Dates.Millisecond(Dates.Hour(1))


function TimingInit()
    ###Construcor. Set the clock and such things
    StartWallTime=now()
    NextWallTime=now()+sec_len ###Wait for a second before showing the clock
    return (StartWallTime,NextWallTime)
end

function get_time_and_unit(Time;Decending=false)
    if Time >= hr_len ###More than an hour
        #println("$Time is bigger than an hour")
        TimeValue=round(Time.value/hr_len.value,digits=1)
        if Decending && Time<(2*hr_len)
            return (TimeValue," hr",Time-hr_len)
        else
            return (TimeValue," hr",hr_len)
        end
    elseif Time >= min_len ###More than a minute
        #println("$Time is bigger than a minute")
        TimeValue=round(Time.value/min_len.value,digits=1)
        if Decending && Time<(2*min_len)
            return (TimeValue," min",Time-min_len)
        else
            return (TimeValue," min",min_len)
        end
    else
        #println("$Time is on the order of seconds")
        return (Int(round(Time.value/sec_len.value))," sec",sec_len)
    end
    return TimeValue, TimeStr, TimeUnit
end

        
function TimingProgress(TimingObject,Nr,NSamples;Message=nothing)
    (StartWallTime,NextTime)=TimingObject
    NowWallTime=now()
    #println("Nr=$Nr NSamples=$NSamples")
    #println("Wait Time: ",WaitTime)
    #println("Time for next write: ",NextTime)
    #println("Clock Time now:      ",NowWallTime)
    if NowWallTime > NextTime
        #println(".................")
        #println("Time to write!")
        TimeWallElapsed=Dates.Millisecond(NowWallTime-StartWallTime)
        #println("Wall Time Elapsed: ",TimeWallElapsed)
        ##Computes samples left
        SamplesLeft = NSamples - Nr
        #println("SamplesLeft=$SamplesLeft")
        ##Compute time left
        TimeWallLeft = Dates.Millisecond(ceil((SamplesLeft/(1.0*Nr))*TimeWallElapsed.value))
        #println("TimeWallLeft=$TimeWallLeft")
        (t_el,dim_el,step_el)=get_time_and_unit(TimeWallElapsed)
        (t_rem,dim_rem,step_rem)=get_time_and_unit(TimeWallLeft,Decending=true)
        (t_tot,dim_tot,step_tot)=get_time_and_unit(TimeWallElapsed+TimeWallLeft)
        WaitTime=maximum([minimum([step_el,step_rem,step_tot]),sec_len])
        #println(Nr," done ", SamplesLeft," to go. Time: ",TimeWallElapsed,
        #        " Time left: ",TimeWallLeft,
        #        " of ",TimeWallElapsed+TimeWallLeft)
        if Message!=nothing
            println(Message)
        end
        println(Nr," done ", SamplesLeft," to go. Time: ",t_el,dim_el,
                ". Time left: ",t_rem,dim_rem,
                " of ",t_tot,dim_tot,".")
        #println("Wait Time: ",WaitTime)
        NextTime=now()+WaitTime
        #println("Time for new next write: ",NextTime)
        return (StartWallTime,NextTime)
    else
        #println("Skip writing")
        return TimingObject
    end
end



TabLevel=TabLevel[1:end-4]
println(TabLevel*"Close ProgressTimer.jl")
