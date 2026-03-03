# Warsztaty-Sztucznej-Inteligencji-2025-2026
Repozytorium zawierające eda i feature selection do późniejszego trenowania modeli w celu przewidzenia aktywności chemicznej związków.

Nie jest to ostateczna wersja, ale jest przeprowadzana tylko na części danych.

Zarówno końcowe cechy jak i ich forma mogą się zmienić w późniejszych fazach pracy.

Repozytorium zawiera:
- Plik "wstepne czyszczenie" - używając web api pobiera dane i wstępnie je filtruje. Wyniki są zapisywane w plikac z rozszerzeniem csv.
- "activity.py", "molecule.py", "target.py" - zmieniają typy niektórych kolumn i redukują ilość cech. Zgodnie z moją wiedzą cechy które zostały odrzucone nie będą przydatne w tym konkretnym zadaniu. Wyniki są zapisywane w plikac z rozszerzeniem csv.
- "dataset.py" - łączy cechy które mają wpływ na aktywność chemiczną w jeden dataset. Prawdopodobnie ostatechne cechy będą się różnić od obecnych.
