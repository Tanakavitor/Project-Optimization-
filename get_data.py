import yfinance as yf

tickers = [
    "AAPL","AMGN","AMZN","AXP","BA","CAT","CRM","CSCO","CVX","DIS",
    "GS","HD","HON","IBM","JNJ","JPM","KO","MCD","MMM","MRK",
    "MSFT","NKE","NVDA","PG","SHW","TRV","UNH","V","VZ","WMT"
]

df = yf.download(
    tickers,
    start="2024-08-01",
    end="2024-12-31",
    interval="1d",
    auto_adjust=True,
    progress=False
)["Close"]

df.to_csv("precos_dowjones_2024H2.csv", index=True)