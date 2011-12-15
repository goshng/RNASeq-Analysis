paste x* | awk '{for(i=t=0;i<NF;) t+=$++i; $0=t}1'
