Реализовано:
	База (15 Баллов)
		-Локальное освещение по модели Фонга
		-Тени
		-Зеркальные отражения
		-Используется >3 материалов
		-Используется 2 или более различных примитивов
		-В сценах 1 и более источников света
		-Рендеринг каждой сцены занимает меньше секунды
		-Разрешение выше 512х512
	Сцены
		Сцена 1
			Примитивы в сцене
				1) Бесконечная плоскость (База)
				2) Сферы (База)
				3) Прямоугольник
			Материалы в сцене
				1) Матовый серый
				2) Полупрозрачный фиолетовый
				3) преломляющая сфера (+1 Балл)
				4) зеркальная сфера
				5) Цветные глянцевые сферы
			Текстуры (+1 Балл)
				1) Плоскость
				2) Шар
				3) Задний фон (не 3d карта)
			Используется многопоточность (+2 Балла)
			Итого: 15 + 1 + 1 + 2 + (Субъективная реалистичность сцены?) + 
				(Непредусмотренный бонус проверяющего?) >= 19

		Сцена 2
			Примитивы в сцене
				1) Сферы (База 1)
				2) Бесконечная плоскость (База)
				3) Прямоугольник
				4) Произвольная 3d модель в виде треугольных мешей (>2000) (+4 Балла) 
				3 различных материала в кольцах сатурна
			Материалы в сцене
				1) Матовый серый
				2) Полупрозрачный фиолетовый
				3) преломляющая сфера (+1 Балл)
				4) зеркальная сфера
			Текстуры в сцене (+1 Балл)
				1) Плоскость
				2) Шар
				3) Задний фон (не 3d карта)
			Используется многопоточность (+2 Балла)
			Итого: 15 + 1 + 1 + 2 + 1 + (Субъективная реалистичность сцены?) + 
				(Непредусмотренный бонус проверяющего?) >= 23

		Сцена 3
			Примитивы в сцене
				1) Сферы (База 1)
			Материалы в сцене
				1) Глянцевый шар
				2) Зеркальный шар
				3) Преломляющий шар (+1 Балл)
			Используется многопоточность (+2 Балла)
			Использована сферическая карта окружения (+1 Балл)
				Использована сферическая панорама комнаты
			Итого 15 + 1 + 2 + 1 + (Субъективная реалистичность сцены?) +
                (Непредусмотренный бонус проверяющего?) >= 19

	Итог по всем сценам
		-База (+15 Баллов)
		-Текстуры (+1 Балл)
		-Использование дополнительных геометрических примитивов (+2 Балла)
		-Использование карт окружения (+1 Балл)
		-Использование многопоточности (+2 Балла)
		-Преломление (+1 Балл)
		-Использование произвольной 3d модель в виде треугольных мешей (>2000) (+4 Балла) 
		-Субъективная реалистичность сцены?
        -Непредусмотренный бонус проверяющего?
	______________________________________________________________________
	15 + 1 + 1 + 2 + 1 + 2 + 4 (Субъективная реалистичность сцены?) + 
	(Непредусмотренный бонус проверяющего?) >= 26 Балла


В ходе выполнения задания были изучены следующие ресурсы
	https://habr.com/ru/post/436790/
	http://ray-tracing.ru/
	https://www.scratchapixel.com/code.php?id=10&origin=/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes&src=0

