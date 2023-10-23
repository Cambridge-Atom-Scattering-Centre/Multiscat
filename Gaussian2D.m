function I = Gaussian2D(x, y, mu, sigma, A)
    I = A*(1/(sigma*sqrt(2*pi)))*exp(-(x - mu(1)).^2/(2*sigma^2)).*exp(-(y - mu(2)).^2/(2*sigma^2));
end