# MolLoom — карта проекта

- Дистрибутив Python называется `molloom`.
- Импортный namespace `polymer_lib` сохранён для совместимости.

## Код

- `polymer_lib/monomer.py` — создание мономеров и исходных 3D-конформеров.
- `polymer_lib/polymerizer.py` — конфигурация, тактичность, выбор конформеров,
  рост цепи, повторы и результат сборки.
- `polymer_lib/poly.py` — сборка RDKit-молекул и Python-проверки геометрии.
- `polymer_lib/_rust_growth.py` и `polymer_lib/rust_growth/` — FFI и опциональное
  Rust-ядро роста цепи.
- `polymer_lib/optimizers/` — оптимизаторы MMFF и OpenFF.
- `polymer_lib/push_off.py` — финальная релаксация цепи через OpenMM.
- `polymer_lib/packer.py`, `_packing_worker.py` и `_rigid_compression.py` —
  упаковка цепей и обработка конфликтов.
- `polymer_lib/api.py` и `runner.py` — JSON API и JSON Lines runner.

## Проверки

- `tests/test_polymerization.py` — рост цепей, кополимеры, геометрия и тактичность.
- `tests/test_openmm_push_off.py` — OpenMM push-off и сериализация.
- `tests/test_packer.py` и `tests/test_rigid_compression.py` — упаковка и сжатие.
- `tests/test_json_api.py` — контракт JSON API.
