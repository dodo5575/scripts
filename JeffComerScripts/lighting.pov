global_settings {
  radiosity {
    pretrace_start 0.2
    pretrace_end   0.2
    count 50

    nearest_count 5
    error_bound 0.5
    recursion_limit 1

    low_error_factor 0.5
    gray_threshold 0.0
    minimum_reuse 0.015
    brightness 0.9

    adc_bailout 0.01/2

  }
}

light_source { 
  <-0.9000, 2.0000, -0.8> 
  color rgb<0.4, 0.4, 0.4> 
  parallel 
  point_at <0.0, 0.0, 0.0> 
}
