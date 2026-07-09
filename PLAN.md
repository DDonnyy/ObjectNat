# ObjectNat 2.0 — план работ

Ветка: `feature/objectnat2.0` → `master`. Целевая версия: **2.0.0** (текущая 1.4.1).

Главная цель: уйти с `networkx` на модель `UrbanGraph` из `iduedu>=2.0.0`, привести
репозиторий и документацию к тому же виду, что и [IduEdu](https://github.com/DDonnyy/IduEdu).

Легенда: `[ ]` не сделано · `[~]` в работе · `[x]` готово · `[?]` нужно решение

---

## 1. Миграция на UrbanGraph

### 1a. Accessibility (изохроны + покрытия) — порт LPRP ✅ написан (2026-07-09)

Новый пакет `objectnat/methods/accessibility/` (зеркало LPRP, импорты из `iduedu`):
- [x] `_utils.py` — geometry-хелперы (voronoi cells, radius/ways/separate, stepped). Self-contained
- [x] `coverage.py` — `get_graph_coverage`, `get_stepped_graph_coverage`
- [x] `isochrones.py` — `get_graph_isochrones`, `get_stepped_graph_isochrones`
- [x] `radius.py` — `get_radius_coverage` (graph-free, перенесён без изменений)
- [x] **GAP 2 fix**: `road_edges_for_ways()` — «ways» на intermodal/walk берёт только `type=="walk"`
- [x] **Решения приняты**: изохроны возвращают только `GeoDataFrame` (без троек ОТ); публичный API — стиль LPRP
- [x] Проверено: файлы компилируются; возвращаемые типы iduedu-функций сверены с использованием в коде
- [x] **Runtime-валидация — 14/14 accessibility-тестов зелёные** на реальном intermodal-графе (iduedu 2.0 + pandas 3.0)

### 1b. Wiring ✅ выполнено и провалидировано на реальном графе (2026-07-09)

- [x] `pyproject.toml`: `iduedu>=2.0.0` в main; поднято `pandas>=3.0.0`, `numpy>=2.4.0`,
      `scipy>=1.17.0`, `geopandas>=1.1.0`; `numba` пришёл транзитивно; `uv lock` (iduedu 2.0.0,
      pandas 3.0.3, numpy 2.4.6, numba 0.65.1). `networkx` пока остаётся — нужен `utils/graph_utils` (§1c)
- [x] `_api.py` / `__init__.py` — экспортируют `get_graph_isochrones`, `get_stepped_graph_isochrones`,
      `get_graph_coverage`, `get_stepped_graph_coverage`, `get_radius_coverage`; старые имена убраны
- [x] Удалены `methods/isochrones/` и `methods/coverage_zones/` (весь networkx-код изохрон/покрытий)
- [x] `tests/conftest.py` — граф на `UrbanGraph`, кэш через `.urbangraph` (write/read), бэкенд `Agg`
- [x] Переписаны `test_isochrones.py`, `test_coverage_zones.py` под новый API (без троек ОТ, `geometry_type`)
- [x] Починен путь к данным в conftest (`../../docs` → `../docs`, поехал после переноса `src/tests`)
- [x] **`import objectnat` под pandas 3.0 — OK**; **14/14 accessibility-тестов зелёные** на реальном
      intermodal-графе (Overpass), PNG перерендерились

- [x] **Хвост починен:** `get_air_resist_ratio` экспортирован из `methods.noise`, импорт в
      `test_noise_simulation.py` поправлен; noise-тесты 5/5 зелёные под pandas 3.0.

### 1c. Остальные утилиты ✅ (2026-07-09)

- [x] `gdf_to_graph` удалён (не нужен); мёртвый кластер `reverse_graph`,
      `remove_weakly_connected_nodes`, `get_closest_nodes_from_gdf` удалён (юзали снесённые модули)
- [x] `math_utils` удалён вместе с мёртвым закомментированным блоком `visibility_analysis` (~230 строк)
- [~] `graph_to_gdf` **оставлен** (публичный, снести не просили) — он последний потребитель `networkx`.
      Решить: убирать ли `graph_to_gdf` + `networkx` целиком (iduedu-`UrbanGraph` отдаёт `edges_gdf`/`nodes_gdf` нативно)

---

## 2. Переезд с poetry на uv ✅ (2026-07-09)

Зеркалим IduEdu: `hatchling` + `uv` + PEP 735 `[dependency-groups]` + `python-semantic-release`.

- [x] `[build-system]` → `hatchling`, версия `dynamic` из `objectnat/_version.py`
- [x] `[tool.poetry]` → `[project]` (PEP 621); `scipy` добавлен явной зависимостью
- [x] Группы `[tool.poetry.group.*]` → `[dependency-groups]`:
      `dev`, `docs`, `lint`, `notebooks`, `release`, `test`, `viz`
- [x] Сгенерирован и коммитится `uv.lock` (208 пакетов); `/poetry.lock` убран из `.gitignore`
- [x] `Makefile`: полностью на `uv` (mirror IduEdu), таргеты `coverage-xml`, `version-next`, `changelog`
- [x] `pylint`: убран `extension-pkg-allow-list = ["networkit"]`
- [x] `scripts/sync_version.py` удалён; версия — единственный источник `objectnat/_version.py`,
      bump через `python-semantic-release`
- [x] CI разбит на 3 воркфлоу 1:1 с IduEdu (`IDUclub/IduEdu`→`IDUclub/ObjectNat`, `main`→`master`):
      `quality.yml` (Tests and Coverage), `release.yml`, `docs.yml`. Старый `ci_pipeline.yml` удалён.
- [x] `quality.yml` — uv, `format-check` → `coverage-xml` → codecov v5 (OIDC); в конце сохранена
      **ObjectNat-специфика**: рендер тест-картинок → force-push в orphan-ветку `assets`
- [x] `release.yml` → `python-semantic-release`, чейнится с `Tests and Coverage` на `master`,
      guard `IDUclub/ObjectNat`, PyPI Trusted Publishing
- [x] `docs.yml` — на PR (paths-filter) строит артефакт, после `Release` деплоит на gh-pages
      (uv + `sphinx-build --keep-going`; `docs/requirements.txt` не нужен — группа `docs`)
- [x] Проверено: `uv build` (1.4.1 из `_version.py`), `import objectnat`, `format-check`,
      `semantic-release --noop version` — всё зелёное

**Осталось на потом (не блокирует):** `CHANGELOG.md`/`CONTRIBUTING.md` (§7), матрица 3.11/3.12 в CI (§4),
markers для network-тестов (§5). Релизы поедут только с ветки `master` при Conventional Commits.

---

## 3. Плоская структура (`src/` → корень) ✅ (2026-07-09)

- [x] `src/objectnat/` → `objectnat/`, `src/tests/` → `tests/` (сделал пользователь)
- [x] `tests/__init__.py` удалён (пакетом быть не должен)
- [x] `pyproject.toml` — hatch `packages = ["objectnat"]` (в §2)
- [x] `Makefile` — `objectnat tests` (в §2); CI — `pytest`/`pylint objectnat` (в §2)
- [x] `docs/conf.py` — `sys.path` `../src` → `..`
- [x] Проверено: ни одного `../src`/`from ...src` не осталось

---

## 4. Починка CI/CD

Найденные баги:

- [x] **Кэш venv никогда не инвалидируется.** Ключ на `poetry.lock` (gitignored) → CI тянул
      протухший venv. Снято в §2: uv-based CI, кэш не нужен.
- [x] Перевести все jobs с `abatilo/actions-poetry` на uv — сделано в §2
- [x] `release.yml`: ручной `tag_version_check` → `python-semantic-release` — сделано в §2
- [ ] **Линтер не гоняется в CI.** Раньше было `pylint ... || true` (не мог упасть); теперь,
      как в IduEdu, CI гоняет только `format-check`, а `pylint` — локально (`make lint`).
      Решить: добавлять ли pylint-гейт в CI и доводить код до прохождения.
- [ ] Прогонять тесты на матрице 3.11 / 3.12 (сейчас один `PY_VER`)

---

## 5. Тесты

Сейчас 1086 строк, `test_utils.py` — 16 строк на два ассерта.

- [x] Переписать фикстуры `conftest.py` под `UrbanGraph` (сделано в §1b: `.urbangraph`-кэш, `Agg`-бэкенд)
- [ ] Замерить базовое покрытие, зафиксировать порог в CI (`--cov-fail-under`)
- [ ] Догнать покрытие:
  - [ ] `utils/geom_utils.py` — сейчас не покрыт напрямую
  - [ ] `utils/graph_utils.py` — после миграции
  - [ ] `provision/provision_exceptions.py` — 4 класса исключений, есть тесты только на 2
  - [ ] `_config.py` — `change_logger_lvl`, `set_enable_tqdm`
- [x] `.gitignore`: добавлены `tests/test_output/` и `tests/test_cache/` (в §7)
- [ ] `[?]` Тесты ходят в Overpass API за графом (`osm_id=1114252`) — вынести в
      `@pytest.mark.network` и кэшировать фикстуру-граф в репозитории?

---

## 6. Документация

Найденные ошибки:

- [ ] **`docs/index.rst` документирует несуществующие функции.**
      `objectnat.get_visibility_accurate` и `objectnat.get_visibilities_from_points`
      исчезли при рефакторинге видимости (`f01e734`, слияние `feature/visibility_refactor`),
      остался только унифицированный `get_visibility`. Ссылки битые.
- [ ] **Quickstart в `docs/index.rst` не запускается:** импортируется
      `get_accessibility_isochrones`, а вызывается `get_accessibility_isochrone_stepped`;
      переменные `point` и `G_intermodal` не определены (граф назван `G`); сломанный отступ.
- [ ] **Markdown внутри reStructuredText.** В `index.rst` блок установки обёрнут в
      ` ``` `, а «Contributions are very welcome» — в `>`; в rst это рендерится буквально.
- [ ] **Битые бейджи.** CI указывает на `IDUclub/ObjectNat`, codecov — на `DDonnyy/ObjectNat`.
      Определиться с каноническим репозиторием.
- [ ] Убрать NetworkX из списка «интегрируется с» в `index.rst` и `README`
- [ ] `README.rst` + `README_RU.rst` → `README.md` (как в IduEdu); решить судьбу русской версии
- [ ] Насытить `docs/methods/examples/` ноутбуками: сейчас 7 штук, все на nx-графах — переписать
- [ ] Добавить страницу «Migration 1.x → 2.0» с таблицей переименований
- [ ] Пересобрать `docs/_build/` не коммитить (уже в `.gitignore`? — **нет, проверить**)

---

## 7. Приведение репозитория в порядок — по аналогии с IduEdu ✅ частично (2026-07-09)

- [x] `AGENTS.md` — полноценный (структура проекта, команды, конвенции, релизы), формат IduEdu
- [x] `CHANGELOG.md` — seed-заголовок, дальше ведёт `python-semantic-release`
- [x] `CONTRIBUTING.md` — адаптирован из IduEdu (IDUclub/ObjectNat, `objectnat/_version.py`, `master`)
- [x] `.editorconfig`, `.gitattributes` (LF + `*.urbangraph binary`) — скопированы из IduEdu
- [x] `.pre-commit-config.yaml` — revs подняты под IduEdu (black 25.9.0, isort 6.1.0). **pylint не добавляю:**
      в IduEdu pre-commit только black+isort (по аналогии)
- [x] Удалён мёртвый закомментированный блок в `visibility_analysis` (в §1c)
- [x] `graphify-out/` **в `.gitignore`** — не коммитим (каждый генерит локально; в IduEdu тоже уйдёт в ignore).
      `AGENTS.md` рекомендует graphify + как активировать (`pip install graphifyy; graphify .`), ссылка на репо
- [x] `.gitignore`: добавлены `tests/test_output/`, `tests/test_cache/`, `graphify-out/`
- [ ] `dev/` — три ноутбука (`iduedutest`, `isochrones_bug`, `noise`), в `.gitignore`. В `docs/examples` или удалить
- [ ] `README.rst`/`README_RU.rst` → `README.md` (§6, как в IduEdu)
- [x] Перегенерить `graphify-out/` под новую структуру (286 узлов / 17 сообществ, отражает accessibility)

---

## 8. Релиз

- [ ] Версия `2.0.0`
- [ ] Migration guide в CHANGELOG: `nx_graph` → `graph`, удалённые `graph_to_gdf`/`gdf_to_graph`,
      удалённый `math_utils`
- [ ] Тег `v2.0.0`, публикация на PyPI

---

## Порядок выполнения

```
Gap #1–#3 в iduedu  ─┐
                     ├─► §1 миграция ─► §5 тесты ─► §6 документация ─► §8 релиз
§2 uv ─► §3 плоская ─┘                    │
                                    §4 CI ─┘
§7 порядок в репо — параллельно, в любой момент
```

Первый шаг, который ничего не блокирует и всё разблокирует: **§2 + §3**
(uv + плоская структура), потому что после них CI начнёт честно падать
и покажет реальное состояние кода перед миграцией.
