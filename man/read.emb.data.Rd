\name{read.emb.data}
\alias{read.emb.data}
\title{Read data from new Spiroff format}

\description{
  Read 2d data from file and return data frame

  \url{http://urchin.spbcas.ru/flyex/apt.jsp?by=name&data=qdata&htmt=mf2iq.htmt&embryo_name=bd3&level=1&order_cl=e.name&button_name=Search}
}

\usage{
read.emb.data(file)
}

\arguments{
  \item{file}{filename}
}

\details{
  Спиров как-то так описывал процедуру
  У них там все сложно и избыточно. По этой опции у меня оно выдает результаты как html стр.
  Насколько я помню / понимаю, я обычно по этой ссылке скачиваю.
  \url{http://urchin.spbcas.ru/flyex/apt.jsp?htmt=mf3_d.htmt&by=name&data=d_load}
  (По опции download)
  и далее quantitative data, ASCII.
  В итоге оно направляет на стр. для скачивания и в конце-концов удается скачать архив с данными: прикладываю.

  ASCII format, where column names are given in a specific form in the first file rows.
  Далее - ссылка, откуда они такие.


  Файл этого формата --- просто числа в столбцах, разделенные пробелами.
  Но почему-то у него нет нормального заголовка, а вместо него шаблонное описание
  столбцов. Эта функция парсит это описание и возвращает датафрейм с корректными именами.
}

\value{
  A data frame.
}

\examples{
xlim <- c(22, 88)
ylim <- c(32, 68)
L <- c(15, 15)

file <- system.file("extdata/data", "ab16.txt", package = "BioSSA")
df <- read.emb.data(file)

bs <- BioSSA(cad ~ AP + DV, data = df, ylim = ylim, xlim = xlim, L = L)
nm <- noise.model(bs, 1:3)

plot(nm)
summary(nm)

# More breif
plot(bs, type = "residuals", groups = 1:3)
}
