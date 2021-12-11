function f=midi2hz(m)
    f= 440 * exp ((m-69) * log(2)/12);
end