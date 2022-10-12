# Chl-binding
This script identify the chl binding sites that have hydrogen bond donors within
a certain (radius) around the (pivot) atom. DN is the list of atoms that are considered
hydrogen bond donor and AN is the list of hydrogen bond acceptor. There is an option to
exclude the the backbone nitrogen in alpha helcies (exclude_alpha_helix) because they
make hydrogen with the oxygens. With the angleT option you may exclude hydrogen
bonds that are on the top of the chl plan, which is defined by the Mg and the N atoms
