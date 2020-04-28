
plot_contour(hours, Zprofile, Tprofile, t_plot, filename) = nothing
get_plot(grid, Q, aux, t) = nothing
export_plots(plots, filename) = nothing

using Requires
@init @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
  using .Plots


  function plot_contour(hours, Zprofile, Tprofile, t_plot, filename)
    contour(hours, Zprofile.*100, Tprofile,
        levels=243.15:323.15, xticks=0:12:t_plot, xlimits=(12,t_plot),
        xlabel="Time of day (hours)", ylabel="Soil depth (cm)", title="Soil temperature (°K)")
    savefig(filename)
  end

  function get_plot(grid, Q, aux, t)
      # TODO:
      # this currently uses some internals: provide a better way to do this
      gridg = get_z(grid)
      Tg = get_data(grid, aux.data)
      return plot(Tg, gridg, ylabel="depth (cm) at t=$(t)", xlabel="T (°K)", yticks=-100:20:0, xlimits=(263.15,303.15), legend=false)
  end

  function export_plots(plots, filename)
      plots isa Tuple || (plots = (plots,))
      plot(plots...)
      savefig(filename)
  end

end
