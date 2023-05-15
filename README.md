W pliku 01 znajduja sie 2 implementacje GLASSO. Pierwsza, `z_zajec_GLASSO`, jest algorytmem który był na zajęciach. Niestety on rozbiegał zamiast zbiegać. Druga, `z_papiera_GLASSO`, jest implementacją z papiera tak jak ją zrozumieliśmy. Skorzystalismy w niej z zbudowanej funkcji na rozwiązywanie LASSO, bo podejrzewaliśmy, że to tam właśnie był problem w `z_zajec_GLASSO`.

Można zawołać funckję `z_papiera_GLASSO(verbose = TRUE)`.

Niestety na danych z p > 50 liczy się dość powoli... Dlatego nie odpalaliśmy nawet na p = 500 z danych 1.
